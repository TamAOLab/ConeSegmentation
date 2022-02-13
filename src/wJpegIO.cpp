#ifdef _WIN32

#define _CRT_SECURE_NO_WARNINGS

#include <sys/stat.h>
#include <locale>
#include <codecvt>
#include <itkVersion.h>
#include <itksys/SystemTools.hxx>
#include <itk_jpeg.h>
#include <setjmp.h>

#include "wJpegIO.h"

// create an error handler for jpeg that
// can longjmp out of the jpeg library
struct itk_jpeg_error_mgr {
	struct jpeg_error_mgr pub;    /* "public" fields */
	jmp_buf setjmp_buffer;        /* for return to caller */
};

extern "C" {
	static void itk_jpeg_error_exit(j_common_ptr cinfo)
	{
		/* cinfo->err really points to a itk_jpeg_error_mgr struct, so coerce pointer
		  */
		itk_jpeg_error_mgr *myerr = reinterpret_cast<itk_jpeg_error_mgr *>(cinfo->err);

		/* Always display the message. */
		/* We could postpone this until after returning, if we chose. */
		(*cinfo->err->output_message)(cinfo);

		jpeg_abort(cinfo);        /* clean up libjpeg state */
		/* Return control to the setjmp point */
		longjmp(myerr->setjmp_buffer, 1);
	}

	static void itk_jpeg_output_message(j_common_ptr)
	{
	}
}

// Wrap setjmp call to avoid warnings about variable clobbering.
static bool wwrapSetjmp(itk_jpeg_error_mgr & jerr)
{
	if (setjmp(jerr.setjmp_buffer))
	{
		return true;
	}
	return false;
}


// simple class to call fopen on construct and
// fclose on destruct
class wJPEGFileWrapper
{
public:
	wJPEGFileWrapper(const char *const fname, const char *const openMode) :m_FilePointer(ITK_NULLPTR)
	{
		std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
		std::wstring wfn = converter.from_bytes(std::string(fname));
		std::wstring wOpenMode = converter.from_bytes(std::string(openMode));

		// std::cout << "wJPEGFileWrapper() : " << fname << std::endl;

		m_FilePointer = _wfopen(wfn.c_str(), wOpenMode.c_str());
	}

	virtual ~wJPEGFileWrapper()
	{
		if (m_FilePointer != ITK_NULLPTR)
		{
			fclose(m_FilePointer);
		}
	}

	FILE *m_FilePointer;
};

bool wJPEGImageIO::CanReadFile(const char *file)
{
	// First check the extension
	std::string filename = file;

	if (filename == "")
	{
		itkDebugMacro(<< "No filename specified.");
		return false;
	}

	bool                   extensionFound = false;
	std::string::size_type JPEGPos = filename.rfind(".jpeg");
	if ((JPEGPos != std::string::npos)
		&& (JPEGPos == filename.length() - 5))
	{
		extensionFound = true;
	}

	JPEGPos = filename.rfind(".JPEG");
	if ((JPEGPos != std::string::npos)
		&& (JPEGPos == filename.length() - 5))
	{
		extensionFound = true;
	}

	JPEGPos = filename.rfind(".jpg");
	if ((JPEGPos != std::string::npos)
		&& (JPEGPos == filename.length() - 4))
	{
		extensionFound = true;
	}

	JPEGPos = filename.rfind(".JPG");
	if ((JPEGPos != std::string::npos)
		&& (JPEGPos == filename.length() - 4))
	{
		extensionFound = true;
	}

	if (!extensionFound)
	{
		itkDebugMacro(<< "The filename extension is not recognized");
		return false;
	}

	// Now check the file header
	wJPEGFileWrapper JPEGfp(file, "rb");
	if (JPEGfp.m_FilePointer == ITK_NULLPTR)
	{
		return false;
	}

	// read the first two bytes
	unsigned char magic[2];
	int           n = static_cast<int>(fread(magic, sizeof(magic), 1, JPEGfp.m_FilePointer));
	if (n != 1)
	{
		return false;
	}

	// check for the magic stuff:
	// 0xFF followed by 0xD8
	if (magic[0] != 0xFF
		|| magic[1] != 0xD8)
	{
		return false;
	}
	// go back to the start of the file
	fseek(JPEGfp.m_FilePointer, 0, SEEK_SET);
	// magic number is ok, try and read the header
	struct itk_jpeg_error_mgr     jerr;
	struct jpeg_decompress_struct cinfo;
	cinfo.err = jpeg_std_error(&jerr.pub);
	// for any jpeg error call itk_jpeg_error_exit
	jerr.pub.error_exit = itk_jpeg_error_exit;
	// for any output message call itk_jpeg_output_message
	jerr.pub.output_message = itk_jpeg_output_message;
	// set the jump point, if there is a jpeg error or warning
	// this will evaluate to true
	if (wwrapSetjmp(jerr))
	{
		// clean up
		jpeg_destroy_decompress(&cinfo);
		// this is not a valid jpeg file
		return false;
	}
	/* Now we can initialize the JPEG decompression object. */
	jpeg_create_decompress(&cinfo);
	/* Step 2: specify data source (eg, a file) */
	jpeg_stdio_src(&cinfo, JPEGfp.m_FilePointer);
	/* Step 3: read file parameters with jpeg_read_header() */
	jpeg_read_header(&cinfo, TRUE);

	// if no errors have occurred yet, then it must be jpeg
	jpeg_destroy_decompress(&cinfo);

	return true;
}

//-----------------------------------------------------------------------------

void wJPEGImageIO::Read(void *buffer)
{
	unsigned int ui;

	// use this class so return will call close
	wJPEGFileWrapper JPEGfp(this->GetFileName(), "rb");
	FILE *          fp = JPEGfp.m_FilePointer;

	if (!fp)
	{
		itkExceptionMacro("Error JPEGImageIO could not open file: "
			<< this->GetFileName()
			<< std::endl
			<< "Reason: "
			<< itksys::SystemTools::GetLastSystemError());
	}

	// create jpeg decompression object and error handler
	struct jpeg_decompress_struct cinfo;
	struct itk_jpeg_error_mgr     jerr;

	cinfo.err = jpeg_std_error(&jerr.pub);
	// for any jpeg error call itk_jpeg_error_exit
	jerr.pub.error_exit = itk_jpeg_error_exit;
	// for any output message call itk_jpeg_output_message
	jerr.pub.output_message = itk_jpeg_output_message;
	if (wwrapSetjmp(jerr))
	{
		// clean up
		jpeg_destroy_decompress(&cinfo);
		itkExceptionMacro("libjpeg could not read file: "
			<< this->GetFileName());
		// this is not a valid jpeg file
	}

	jpeg_create_decompress(&cinfo);

	// set the source file
	jpeg_stdio_src(&cinfo, fp);

	// read the header
	jpeg_read_header(&cinfo, TRUE);

	// prepare to read the bulk data
	jpeg_start_decompress(&cinfo);

	SizeValueType rowbytes = cinfo.output_components * cinfo.output_width;
	JSAMPLE *     tempImage = static_cast<JSAMPLE *>(buffer);

	JSAMPROW *row_pointers = new JSAMPROW[cinfo.output_height];
	for (ui = 0; ui < cinfo.output_height; ++ui)
	{
		row_pointers[ui] = tempImage + rowbytes * ui;
	}

	// read the bulk data
	unsigned int remainingRows;
	while (cinfo.output_scanline < cinfo.output_height)
	{
		remainingRows = cinfo.output_height - cinfo.output_scanline;
		jpeg_read_scanlines(&cinfo, &row_pointers[cinfo.output_scanline],
			remainingRows);
	}

	// finish the decompression step
	jpeg_finish_decompress(&cinfo);

	// destroy the decompression object
	jpeg_destroy_decompress(&cinfo);

	delete[] row_pointers;
}

void wJPEGImageIO::ReadImageInformation()
{
	m_Spacing[0] = 1.0;  // We'll look for JPEG pixel size information later,
	m_Spacing[1] = 1.0;  // but set the defaults now

	m_Origin[0] = 0.0;
	m_Origin[1] = 0.0;

	// use this class so return will call close
	wJPEGFileWrapper JPEGfp(m_FileName.c_str(), "rb");
	FILE *          fp = JPEGfp.m_FilePointer;
	if (!fp)
	{
		itkExceptionMacro("Error wJPEGImageIO could not open file: "
			<< this->GetFileName()
			<< std::endl
			<< "Reason: "
			<< itksys::SystemTools::GetLastSystemError());
	}

	// create jpeg decompression object and error handler
	struct jpeg_decompress_struct cinfo;
	struct itk_jpeg_error_mgr     jerr;

	cinfo.err = jpeg_std_error(&jerr.pub);
	jerr.pub.error_exit = itk_jpeg_error_exit;
	if (setjmp(jerr.setjmp_buffer))
	{
		// clean up
		jpeg_destroy_decompress(&cinfo);
		// this is not a valid jpeg file
		itkExceptionMacro("Error wJPEGImageIO could not open file: "
			<< this->GetFileName());
	}
	jpeg_create_decompress(&cinfo);

	// set the source file
	jpeg_stdio_src(&cinfo, fp);

	// read the header
	jpeg_read_header(&cinfo, TRUE);

	// force the output image size to be calculated (we could have used
	// cinfo.image_height etc. but that would preclude using libjpeg's
	// ability to scale an image on input).
	jpeg_calc_output_dimensions(&cinfo);

	// pull out the width/height
	this->SetNumberOfDimensions(2);
	m_Dimensions[0] = cinfo.output_width;
	m_Dimensions[1] = cinfo.output_height;

	this->SetNumberOfComponents(cinfo.output_components);

	switch (this->GetNumberOfComponents())
	{
	case 1:
		m_PixelType = SCALAR;
		break;
	case 2:
		m_PixelType = VECTOR;
		break;
	case 3:
		m_PixelType = RGB;
		break;
	case 4:
		m_PixelType = RGBA;
		break;
	}

	// If we have some spacing information we use it
	if (cinfo.density_unit > 0
		&& cinfo.X_density > 0
		&& cinfo.Y_density > 0
		)
	{
		if (cinfo.density_unit == 1) // inches
		{
			m_Spacing[0] = 25.4 / cinfo.X_density;
			m_Spacing[1] = 25.4 / cinfo.Y_density;
		}
		else if (cinfo.density_unit == 2) // cm
		{
			m_Spacing[0] = 10.0 / cinfo.X_density;
			m_Spacing[1] = 10.0 / cinfo.Y_density;
		}
	}

	// close the file
	jpeg_destroy_decompress(&cinfo);
}

void wJPEGImageIO::Write(const void *buffer)
{
	// the IORegion is not required to be set so we must use GetNumberOfDimensions
	if (this->GetNumberOfDimensions() != 2)
	{
		itkExceptionMacro(<< "JPEG Writer can only write 2-dimensional images");
	}

	if (this->GetComponentType() != UCHAR
		&& this->GetComponentType() != UINT)
	{
		itkExceptionMacro(<< "JPEG supports unsigned char/int only");
	}

	this->WriteSlice(m_FileName, buffer);
}

void wJPEGImageIO::WriteSlice(std::string & fileName, const void *buffer)
{
	// use this class so return will call close
	wJPEGFileWrapper JPEGfp(fileName.c_str(), "wb");
	FILE *          fp = JPEGfp.m_FilePointer;

	if (!fp)
	{
		itkExceptionMacro("Unable to open file "
			<< fileName
			<< " for writing."
			<< std::endl
			<< "Reason: "
			<< itksys::SystemTools::GetLastSystemError());
	}

	// Call the correct templated function for the output

	// overriding jpeg_error_mgr so we don't exit when an error happens
	// Create the jpeg compression object and error handler
	//struct jpeg_compress_struct cinfo;
	//struct itk_jpeg_error_mgr jerr;

	struct itk_jpeg_error_mgr   jerr;
	struct jpeg_compress_struct cinfo;
	cinfo.err = jpeg_std_error(&jerr.pub);
	// set the jump point, if there is a jpeg error or warning
	// this will evaluate to true
	if (wwrapSetjmp(jerr))
	{
		jpeg_destroy_compress(&cinfo);
		itkExceptionMacro(<< "JPEG : Out of disk space");
	}

	jpeg_create_compress(&cinfo);

	// set the destination file
	//struct jpeg_destination_mgr compressionDestination;
	jpeg_stdio_dest(&cinfo, fp);

	// set the information about image
	const SizeValueType width = m_Dimensions[0];
	const SizeValueType height = m_Dimensions[1];

	// The JPEG standard only supports images up to 64K*64K due to 16-bit fields
	// in SOF markers.
	cinfo.image_width = width;
	cinfo.image_height = height;
	if (cinfo.image_width > 65536 || cinfo.image_height > 65536)
	{
		itkExceptionMacro(<< "JPEG : Image is too large for JPEG");
	}

	cinfo.input_components = this->GetNumberOfComponents();
	const unsigned int numComp = this->GetNumberOfComponents();

	// Maximum number of components (color channels) allowed in JPEG image.
	// JPEG spec set this to 255. However ijg default it to 10.
	if (cinfo.input_components > 255)
	{
		itkExceptionMacro(<< "JPEG : Too many components for JPEG");
	}
	if (cinfo.input_components > MAX_COMPONENTS)
	{
		itkExceptionMacro(<< "JPEG : Too many components for IJG. Recompile IJG.");
	}

	switch (cinfo.input_components)
	{
	case 1:
		cinfo.in_color_space = JCS_GRAYSCALE;
		break;
	case 3:
		cinfo.in_color_space = JCS_RGB;
		break;
	default:
		cinfo.in_color_space = JCS_UNKNOWN;
		break;
	}

	// set the compression parameters
	jpeg_set_defaults(&cinfo);         // start with reasonable defaults
	jpeg_set_quality(&cinfo, m_Quality, TRUE);
	if (m_Progressive)
	{
		jpeg_simple_progression(&cinfo);
	}

	if (m_Spacing[0] > 0 && m_Spacing[1] > 0)
	{
		// store the spacing information as pixels per inch or cm, depending on which option
		// retains as much precision as possible
		std::vector< UINT16 > densityPerInch(2);
		densityPerInch[0] = static_cast<UINT16>(25.4 / m_Spacing[0] + 0.5);
		densityPerInch[1] = static_cast<UINT16>(25.4 / m_Spacing[1] + 0.5);

		std::vector< UINT16 > densityPerCm(2);
		densityPerCm[0] = static_cast<UINT16>(10.0 / m_Spacing[0] + 0.5);
		densityPerCm[1] = static_cast<UINT16>(10.0 / m_Spacing[1] + 0.5);

		if (std::abs(25.4 / m_Spacing[0] - densityPerInch[0]) + std::abs(25.4 / m_Spacing[1] - densityPerInch[1])
			<= std::abs(10.0 / m_Spacing[0] - densityPerCm[0]) + std::abs(10.0 / m_Spacing[1] - densityPerCm[1]))
		{
			cinfo.density_unit = 1;
			cinfo.X_density = densityPerInch[0];
			cinfo.Y_density = densityPerInch[1];
		}
		else
		{
			cinfo.density_unit = 0;
			cinfo.X_density = densityPerCm[0];
			cinfo.Y_density = densityPerCm[1];
		}
	}

	// start compression
	jpeg_start_compress(&cinfo, TRUE);

	volatile const JSAMPLE *outPtr = ((const JSAMPLE *)buffer);

	// write the data. in jpeg, the first row is the top row of the image
	JSAMPROW *row_pointers = new JSAMPROW[height];
	const int rowInc = numComp * width;
	for (unsigned int ui = 0; ui < height; ui++)
	{
		row_pointers[ui] = const_cast<JSAMPROW>(outPtr);
		outPtr = const_cast<JSAMPLE *>(outPtr) + rowInc;
	}
	jpeg_write_scanlines(&cinfo, row_pointers, height);

	if (fflush(fp) == EOF)
	{
		itkExceptionMacro(<< "JPEG : Out of disk space");
	}

	// finish the compression
	jpeg_finish_compress(&cinfo);

	// clean up and close the file
	delete[] row_pointers;
	jpeg_destroy_compress(&cinfo);
}



wJPEGImageIOFactory::wJPEGImageIOFactory()
{
	this->RegisterOverride("itkImageIOBase",
		"wJpegIO",
		"JPEG Image IO (Unicode)",
		1,
		itk::CreateObjectFunction< wJPEGImageIO >::New());

	// std::cout << "wJPEGImageIOFactory()" << std::endl;
}

wJPEGImageIOFactory::~wJPEGImageIOFactory()
{}

const char *
wJPEGImageIOFactory::GetITKSourceVersion(void) const
{
	return ITK_SOURCE_VERSION;
}

const char *
wJPEGImageIOFactory::GetDescription(void) const
{
	return "JPEG ImageIO Factory (Unicode), allows the loading of JPEG images into insight";
}

// Automatically register wJPEGImageIOFactory with ITK on startup,
// making wJPEGImageIO selectable for ImageFileReader.
static class wJpegSelfRegistration
{
public:
	wJpegSelfRegistration()
	{
		// Find and un-register the original factory
		std::string oldname = "JPEGImageIOFactory";
		std::list<itk::ObjectFactoryBase *> factories = itk::ObjectFactoryBase::GetRegisteredFactories();
		for (itk::ObjectFactoryBase *f : factories) {
			if (oldname == std::string(f->GetNameOfClass())) {
				// std::cout << "Unregister " << f->GetNameOfClass() << std::endl;
				itk::ObjectFactoryBase::UnRegisterFactory(f);
				break;
			}
		}
		// Now register the unicode version
		wJPEGImageIOFactory::RegisterOneFactory();
	}
} jpegSelfRegistration;


#endif
