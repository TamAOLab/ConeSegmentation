#ifdef _WIN32

#define _CRT_SECURE_NO_WARNINGS

#include <sys/stat.h>
#include <locale>
#include <codecvt>
#include <itkVersion.h>
#include <itk_png.h>
#include <itksys/SystemTools.hxx>

#include "wPngIO.h"

extern "C"
{
#include <setjmp.h>
	/* The PNG library does not expect the error function to return.
	   Therefore we must use this ugly longjmp call.  */
	static void wPNGWriteErrorFunction(png_structp png_ptr,
		png_const_charp itkNotUsed(error_msg))
	{
		longjmp(png_jmpbuf(png_ptr), 1);
	}
}

extern "C"
{
	static void wPNGWriteWarningFunction(png_structp itkNotUsed(png_ptr),
		png_const_charp itkNotUsed(warning_msg))
	{}
}

static bool wwrapSetjmp(png_structp & png_ptr)
{
	if (setjmp(png_jmpbuf(png_ptr)))
	{
		return 1;
	}
	return 0;
}

// simple class to call fopen on construct and
// fclose on destruct
class wPNGFileWrapper
{
public:
	wPNGFileWrapper(const char *const fname, const char *const openMode) :m_FilePointer(ITK_NULLPTR)
	{
		std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
		std::wstring wfn = converter.from_bytes(std::string(fname));
		std::wstring wOpenMode = converter.from_bytes(std::string(openMode));

		// std::cout << "wPNGFileWrapper() : " << fname << std::endl;

		m_FilePointer = _wfopen(wfn.c_str(), wOpenMode.c_str());
	}

	virtual ~wPNGFileWrapper()
	{
		if (m_FilePointer)
		{
			fclose(m_FilePointer);
		}
	}

	FILE *m_FilePointer;
};

bool wPNGImageIO::CanReadFile(const char *file)
{
	// First check the filename
	std::string filename = file;

	if (filename == "")
	{
		itkDebugMacro(<< "No filename specified.");
		return false;
	}

	// std::cout << "wPNGImageIO::CanReadFile() : " << filename << std::endl;

	// Now check the file header
	wPNGFileWrapper pngfp(file, "rb");
	if (pngfp.m_FilePointer == ITK_NULLPTR)
	{
		return false;
	}
	unsigned char header[8];
	size_t temp = fread(header, 1, 8, pngfp.m_FilePointer);
	if (temp != 8)
	{
		return false;
	}
	bool is_png = !png_sig_cmp(header, 0, 8);
	if (!is_png)
	{
		return false;
	}
	png_structp png_ptr = png_create_read_struct
	(PNG_LIBPNG_VER_STRING, (png_voidp)ITK_NULLPTR,
		ITK_NULLPTR, ITK_NULLPTR);
	if (!png_ptr)
	{
		return false;
	}

	png_infop info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
	{
		png_destroy_read_struct(&png_ptr,
			(png_infopp)ITK_NULLPTR, (png_infopp)ITK_NULLPTR);
		return false;
	}

	png_infop end_info = png_create_info_struct(png_ptr);
	if (!end_info)
	{
		png_destroy_read_struct(&png_ptr, &info_ptr,
			(png_infopp)ITK_NULLPTR);
		return false;
	}
	png_destroy_read_struct(&png_ptr, &info_ptr,
		&end_info);

	return true;
}

void wPNGImageIO::Read(void *buffer)
{
	// std::cout << "wPNGImageIO::Read()" << std::endl;

	itkDebugMacro("Read: file dimensions = " << this->GetNumberOfDimensions());
	// use this class so return will call close
	wPNGFileWrapper pngfp(this->GetFileName(), "rb");
	FILE *         fp = pngfp.m_FilePointer;
	if (!fp)
	{
		itkExceptionMacro("wPNGImageIO could not open file: "
			<< this->GetFileName() << " for reading."
			<< std::endl
			<< "Reason: "
			<< itksys::SystemTools::GetLastSystemError());
	}
	unsigned char header[8];
	size_t temp = fread(header, 1, 8, fp);
	if (temp != 8)
	{
		itkExceptionMacro("wPNGImageIO failed to read header for file: "
			<< this->GetFileName() << std::endl
			<< "Reason: fread read only " << temp
			<< " instead of 8");
	}

	bool is_png = !png_sig_cmp(header, 0, 8);
	if (!is_png)
	{
		itkExceptionMacro("File is not png type: " << this->GetFileName());
	}
	png_structp png_ptr = png_create_read_struct
	(PNG_LIBPNG_VER_STRING, (png_voidp)ITK_NULLPTR,
		ITK_NULLPTR, ITK_NULLPTR);
	if (!png_ptr)
	{
		itkExceptionMacro("File is not png type" << this->GetFileName());
	}

	png_infop info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
	{
		png_destroy_read_struct(&png_ptr,
			(png_infopp)ITK_NULLPTR, (png_infopp)ITK_NULLPTR);
		itkExceptionMacro("File is not png type " << this->GetFileName());
	}

	png_infop end_info = png_create_info_struct(png_ptr);
	if (!end_info)
	{
		png_destroy_read_struct(&png_ptr, &info_ptr,
			(png_infopp)ITK_NULLPTR);
		itkExceptionMacro("File is not png type " << this->GetFileName());
	}

	//  VS 7.1 has problems with setjmp/longjmp in C++ code
#if !defined( MSC_VER ) || _MSC_VER != 1310
	if (wwrapSetjmp(png_ptr))
	{
		png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
		itkExceptionMacro("File is not png type " << this->GetFileName());
	}
#endif

	png_init_io(png_ptr, fp);
	png_set_sig_bytes(png_ptr, 8);

	png_read_info(png_ptr, info_ptr);

	png_uint_32 width, height;
	int         bitDepth, colorType, interlaceType;
	int         compression_type, filter_method;
	png_get_IHDR(png_ptr, info_ptr,
		&width, &height,
		&bitDepth, &colorType, &interlaceType,
		&compression_type, &filter_method);

	if (colorType == PNG_COLOR_TYPE_PALETTE)
	{
		if (this->GetExpandRGBPalette())
		{// convert palette to RGB
			png_set_palette_to_rgb(png_ptr);
		}
		else
		{
			// unpack the pixels
			png_set_packing(png_ptr);
		}
	}

	// minimum of a byte per pixel
	if (colorType == PNG_COLOR_TYPE_GRAY && bitDepth < 8)
	{
#if (PNG_LIBPNG_VER_MAJOR < 2 && PNG_LIBPNG_VER_MINOR < 4)
		png_set_gray_1_2_4_to_8(png_ptr);
#else
		png_set_expand_gray_1_2_4_to_8(png_ptr);
#endif
	}

	// add alpha if any alpha found
	if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
	{
		png_set_tRNS_to_alpha(png_ptr);
	}

	if (bitDepth > 8)
	{
#ifndef ITK_WORDS_BIGENDIAN
		png_set_swap(png_ptr);
#endif
	}

#if (PNG_LIBPNG_VER_MAJOR < 2 && PNG_LIBPNG_VER_MINOR < 5)
	if (info_ptr->valid & PNG_INFO_sBIT)
	{
		png_set_shift(png_ptr, &(info_ptr->sig_bit));
	}
#else
	if (png_get_valid(png_ptr, info_ptr, PNG_INFO_sBIT))
	{
		png_color_8p bits;
		png_get_sBIT(png_ptr, info_ptr, &bits);
		png_set_shift(png_ptr, bits);
	}
#endif
	// have libpng handle interlacing
	//int number_of_passes = png_set_interlace_handling(png_ptr);
	// update the info now that we have defined the filters
	png_read_update_info(png_ptr, info_ptr);

	SizeValueType  rowbytes = static_cast<SizeValueType>(png_get_rowbytes(png_ptr, info_ptr));
	unsigned char *tempImage = static_cast<unsigned char *>(buffer);
	png_bytep *    row_pointers = new png_bytep[height];
	for (unsigned int ui = 0; ui < height; ++ui)
	{
		row_pointers[ui] = tempImage + rowbytes * ui;
	}
	png_read_image(png_ptr, row_pointers);
	delete[] row_pointers;
	// close the file
	png_read_end(png_ptr, ITK_NULLPTR);
	png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);

	// std::cout << "wPNGImageIO::Read() -- done!" << std::endl;
}

void wPNGImageIO::ReadImageInformation()
{
	m_Spacing[0] = 1.0;  // We'll look for PNG pixel size information later,
	m_Spacing[1] = 1.0;  // but set the defaults now

	m_Origin[0] = 0.0;
	m_Origin[1] = 0.0;

	// use this class so return will call close
	wPNGFileWrapper pngfp(m_FileName.c_str(), "rb");
	FILE *         fp = pngfp.m_FilePointer;
	if (!fp)
	{
		return;
	}
	unsigned char header[8];
	size_t temp = fread(header, 1, 8, fp);
	if (temp != 8)
	{
		itkExceptionMacro("wPNGImageIO failed to read header for file: "
			<< this->GetFileName() << std::endl
			<< "Reason: fread read only " << temp
			<< " instead of 8");
	}

	bool is_png = !png_sig_cmp(header, 0, 8);
	if (!is_png)
	{
		return;
	}
	png_structp png_ptr = png_create_read_struct
	(PNG_LIBPNG_VER_STRING, (png_voidp)ITK_NULLPTR,
		ITK_NULLPTR, ITK_NULLPTR);
	if (!png_ptr)
	{
		return;
	}

	png_infop info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
	{
		png_destroy_read_struct(&png_ptr,
			(png_infopp)ITK_NULLPTR, (png_infopp)ITK_NULLPTR);
		return;
	}

	png_infop end_info = png_create_info_struct(png_ptr);
	if (!end_info)
	{
		png_destroy_read_struct(&png_ptr, &info_ptr,
			(png_infopp)ITK_NULLPTR);
		return;
	}

	png_init_io(png_ptr, fp);
	png_set_sig_bytes(png_ptr, 8);

	png_read_info(png_ptr, info_ptr);

	png_uint_32 width, height;
	int         bitDepth, colorType, interlaceType;
	int         compression_type, filter_method;
	png_get_IHDR(png_ptr, info_ptr,
		&width, &height,
		&bitDepth, &colorType, &interlaceType,
		&compression_type, &filter_method);

	m_IsReadAsScalarPlusPalette = false;
	if (colorType == PNG_COLOR_TYPE_PALETTE)
	{
		if (m_ExpandRGBPalette)
		{// convert palettes to RGB
			png_set_palette_to_rgb(png_ptr);
		}
		else
		{
			// Unpack the pixels
			png_set_packing(png_ptr);

			m_IsReadAsScalarPlusPalette = true;

			png_colorp palette;
			int num_entry;
			png_get_PLTE(png_ptr, info_ptr, &palette, &num_entry);

			if (num_entry < 0) num_entry = 0;
			size_t num_entryI(static_cast<size_t>(num_entry));

			m_ColorPalette.resize(num_entryI);
			for (size_t c = 0; c < num_entryI; ++c)
			{
				RGBPixelType p;
				p.SetRed(palette[c].red);
				p.SetGreen(palette[c].green);
				p.SetBlue(palette[c].blue);
				m_ColorPalette[c] = p;
			}

		}
	}

	// minimum of a byte per pixel
	if (colorType == PNG_COLOR_TYPE_GRAY && bitDepth < 8)
	{
#if ( PNG_LIBPNG_VER_MAJOR < 2 && PNG_LIBPNG_VER_MINOR < 4 )
		png_set_gray_1_2_4_to_8(png_ptr);
#else
		png_set_expand_gray_1_2_4_to_8(png_ptr);
#endif
	}

	// add alpha if any alpha found
	if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
	{
		png_set_tRNS_to_alpha(png_ptr);
	}

	// update the info now that we have defined the filters
	png_read_update_info(png_ptr, info_ptr);
	this->SetNumberOfDimensions(2);
	m_Dimensions[0] = width;
	m_Dimensions[1] = height;
	if (bitDepth <= 8)
	{
		m_PixelType = SCALAR;
		m_ComponentType = UCHAR;
	}
	else
	{
		m_PixelType = SCALAR;
		m_ComponentType = USHORT;
	}
	this->SetNumberOfComponents(png_get_channels(png_ptr, info_ptr));

	if (this->GetNumberOfComponents() == 3)
	{
		m_PixelType = RGB;
	}
	else if (this->GetNumberOfComponents() == 4)
	{
		m_PixelType = RGBA;
	}

	// see if the PNG file stored spacing information,
	double px_width = 1.0;
	double px_height = 1.0;

#if defined(PNG_sCAL_SUPPORTED) && defined(PNG_FLOATING_POINT_SUPPORTED)
	int    units = PNG_SCALE_UNKNOWN;

	if (PNG_INFO_sCAL == png_get_sCAL(png_ptr, info_ptr, &units, &px_width, &px_height) &&
		units == PNG_SCALE_UNKNOWN &&
		(px_width != 1.0 || px_height != 1.0))
	{
		// Only libpng <1.5 can read sCAL with SCALE_UNKNOWN, warn this is
		// not going to be compatible with newer libpngs
		itkWarningMacro("PNG sCAL SCALE_UNKNOWN detected with non-unit spacing. This is no longer supported by libpng. Re-saving this file is recommended.");
	}
#endif

	m_Spacing[0] = px_width;
	m_Spacing[1] = px_height;

	// clean up
	png_destroy_read_struct(&png_ptr, &info_ptr,
		&end_info);
}

void wPNGImageIO::Write(const void *buffer)
{
	this->WriteSlice(m_FileName, buffer);
}

void wPNGImageIO::WriteSlice(const std::string & fileName, const void *buffer)
{
	// use this class so return will call close
	wPNGFileWrapper pngfp(fileName.c_str(), "wb");
	FILE *         fp = pngfp.m_FilePointer;

	if (!fp)
	{
		// IMPORTANT: The itkExceptionMacro() cannot be used here due to a bug in
		// Visual
		//            Studio 7.1 in release mode. That compiler will corrupt the
		// RTTI type
		//            of the Exception and prevent the catch() from recognizing it.
		//            For details, see Bug #1872 in the bugtracker.

		::itk::ExceptionObject excp(__FILE__, __LINE__, "Problem while opening the file.", ITK_LOCATION);
		throw excp;
	}

	volatile int bitDepth;
	switch (this->GetComponentType())
	{
	case UCHAR:
		bitDepth = 8;
		break;

	case USHORT:
		bitDepth = 16;
		break;

	default:
	{
		// IMPORTANT: The itkExceptionMacro() cannot be used here due to a bug in
		// Visual
		//            Studio 7.1 in release mode. That compiler will corrupt the
		// RTTI type
		//            of the Exception and prevent the catch() from recognizing
		// it.
		//            For details, see Bug #1872 in the bugtracker.
		::itk::ExceptionObject excp(__FILE__, __LINE__, "PNG supports unsigned char and unsigned short", ITK_LOCATION);
		throw excp;
	}
	}

	png_structp png_ptr = png_create_write_struct
	(PNG_LIBPNG_VER_STRING, (png_voidp)ITK_NULLPTR, ITK_NULLPTR, ITK_NULLPTR);
	if (!png_ptr)
	{
		itkExceptionMacro(<< "Unable to write PNG file! png_create_write_struct failed.");
	}

	png_infop info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
	{
		png_destroy_write_struct(&png_ptr,
			(png_infopp)ITK_NULLPTR);
		itkExceptionMacro(<< "Unable to write PNG file!. png_create_info_struct failed.");
	}

	png_init_io(png_ptr, fp);

	//  VS 7.1 has problems with setjmp/longjmp in C++ code
#if !defined( _MSC_VER ) || _MSC_VER != 1310
	png_set_error_fn(png_ptr, png_ptr,
		wPNGWriteErrorFunction, wPNGWriteWarningFunction);
	if (wwrapSetjmp(png_ptr))
	{
		fclose(fp);
		itkExceptionMacro("Error while writing Slice to file: "
			<< this->GetFileName()
			<< std::endl
			<< "Reason: "
			<< itksys::SystemTools::GetLastSystemError());
	}
#endif

	int          colorType;
	unsigned int numComp = this->GetNumberOfComponents();
	switch (numComp)
	{
	case 1:
		colorType = PNG_COLOR_TYPE_GRAY;
		break;
	case 2:
		colorType = PNG_COLOR_TYPE_GRAY_ALPHA;
		break;
	case 3:
		colorType = PNG_COLOR_TYPE_RGB;
		break;
	default:
		colorType = PNG_COLOR_TYPE_RGB_ALPHA;
		break;
	}

	png_uint_32 width, height;
	double      rowSpacing, colSpacing;
	width = this->GetDimensions(0);
	colSpacing = m_Spacing[0];

	if (m_NumberOfDimensions > 1)
	{
		height = this->GetDimensions(1);
		rowSpacing = m_Spacing[1];
	}
	else
	{
		height = 1;
		rowSpacing = 1;
	}

	png_set_IHDR(png_ptr, info_ptr, width, height,
		bitDepth, colorType, PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);
	// interlaceType - PNG_INTERLACE_NONE or
	//                 PNG_INTERLACE_ADAM7

	if (m_UseCompression)
	{
		// Set the image compression level.
		png_set_compression_level(png_ptr, m_CompressionLevel);
	}

	// write out the spacing information:
	//      set the unit_type to unknown.  if we add units to ITK, we should
	//          convert pixel size to meters and store units as meters (png
	//          has three set of units: meters, radians, and unknown).
#if defined(PNG_sCAL_SUPPORTED) && defined(PNG_FLOATING_POINT_SUPPORTED)
	png_set_sCAL(png_ptr, info_ptr, PNG_SCALE_METER, colSpacing,
		rowSpacing);
#endif

	png_write_info(png_ptr, info_ptr);
	// default is big endian
	if (bitDepth > 8)
	{
#ifndef ITK_WORDS_BIGENDIAN
		png_set_swap(png_ptr);
#endif
	}
	png_byte **row_pointers = new png_byte *[height];

	{
		const int        rowInc = width * numComp * bitDepth / 8;
		volatile const unsigned char *outPtr = ((const unsigned char *)buffer);
		for (unsigned int ui = 0; ui < height; ui++)
		{
			row_pointers[ui] = const_cast<png_byte *>(outPtr);
			outPtr = const_cast<unsigned char *>(outPtr) + rowInc;
		}
	}
	png_write_image(png_ptr, row_pointers);
	png_write_end(png_ptr, info_ptr);

	delete[] row_pointers;
	png_destroy_write_struct(&png_ptr, &info_ptr);
}


wPNGImageIOFactory::wPNGImageIOFactory()
{
	this->RegisterOverride("itkImageIOBase",
		"wPngIO",
		"PNG Image IO (Unicode)",
		1,
		itk::CreateObjectFunction< wPNGImageIO >::New());
}

wPNGImageIOFactory::~wPNGImageIOFactory()
{}

const char *
wPNGImageIOFactory::GetITKSourceVersion(void) const
{
	return ITK_SOURCE_VERSION;
}

const char *
wPNGImageIOFactory::GetDescription(void) const
{
	return "PNG ImageIO Factory (Unicode), allows the loading of PNG images into insight";
}

// Automatically register wPNGImageIOFactory with ITK on startup,
// making wPNGImageIO selectable for ImageFileReader.
static class wPngSelfRegistration
{
public:
	wPngSelfRegistration()
	{
		// Find and un-register the original factory
		std::string oldname = "PNGImageIOFactory";
		std::list<itk::ObjectFactoryBase *> factories = itk::ObjectFactoryBase::GetRegisteredFactories();
		for (itk::ObjectFactoryBase *f : factories) {
			if (oldname == f->GetNameOfClass()) {
				// std::cout << "Unregister " << f->GetNameOfClass() << std::endl;
				itk::ObjectFactoryBase::UnRegisterFactory(f);
				break;
			}
		}
		// Now register the unicode version
		wPNGImageIOFactory::RegisterOneFactory();
	}
} pngSelfRegistration;

#endif

