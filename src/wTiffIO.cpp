#ifdef _WIN32

#include <sys/stat.h>
#include <locale>
#include <codecvt>
#include <itkVersion.h>
#include <itksys/SystemTools.hxx>
#include <itkMetaDataObject.h>
#include "itkTIFFReaderInternal.h"

#include "wTiffIO.h"

static int WTiffReaderInternalOpen(itk::TIFFReaderInternal *_this, const char *filename)
{
	std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
	std::string fn(filename);
	std::wstring wfn = converter.from_bytes(fn);

	_this->Clean();
	struct _stat64i32 fs;
	if (_wstat64i32(wfn.c_str(), &fs))
	{
#if defined(_WIN32) && ! defined(__MINGW32_VERSION)
		struct _stat64 fs64;
		if (_wstat64(wfn.c_str(), &fs64))
		{
			// Both stat() and _stat64() return != 0
			return 0;
		}
#else
		return 0;
#endif
	}

	_this->m_Image = TIFFOpenW(wfn.c_str(), "r");
	if (!_this->m_Image)
	{
		_this->Clean();
		return 0;
	}
	if (!_this->Initialize())
	{
		_this->Clean();
		return 0;
	}

	_this->m_IsOpen = true;
	return 1;
}

bool wTIFFImageIO::CanReadFile(const char *file)
{
	// First check the filename
	std::string filename = file;

	if (filename == "")
	{
		itkDebugMacro(<< "No filename specified.");
		return false;
	}
	
	// std::cout << "wTIFFImageIO::CanReadFile() : " << filename << std::endl;

	// Now check if this is a valid TIFF image
	TIFFErrorHandler save = TIFFSetErrorHandler(ITK_NULLPTR);
	int res = WTiffReaderInternalOpen(m_InternalImage, file);
	if (res)
	{
		TIFFSetErrorHandler(save);
		return true;
	}
	m_InternalImage->Clean();
	TIFFSetErrorHandler(save);
	return false;
}

void wTIFFImageIO::Write(const void *buffer)
{
	if (m_NumberOfDimensions == 2 || m_NumberOfDimensions == 3)
	{
		this->InternalWrite(buffer);
	}
	else
	{
		itkExceptionMacro(<< "TIFF Writer can only write 2-d or 3-d images");
	}
}

void wTIFFImageIO::InternalWrite(const void *buffer)
{
	const char *outPtr = (const char *)buffer;

	unsigned int page, pages = 1;

	const SizeValueType width = m_Dimensions[0];
	const SizeValueType height = m_Dimensions[1];
	if (m_NumberOfDimensions == 3)
	{
		pages = m_Dimensions[2];
	}

	int    scomponents = this->GetNumberOfComponents();
	float  resolution_x = static_cast<float>(m_Spacing[0] != 0.0 ? 25.4 / m_Spacing[0] : 0.0);
	float  resolution_y = static_cast<float>(m_Spacing[1] != 0.0 ? 25.4 / m_Spacing[1] : 0.0);
	// rowsperstrip is set to a default value but modified based on the tif scanlinesize before
	// passing it into the TIFFSetField (see below).
	uint32 rowsperstrip = (uint32)-1;
	int    bps;

	switch (this->GetComponentType())
	{
	case UCHAR:
		bps = 8;
		break;
	case CHAR:
		bps = 8;
		break;
	case USHORT:
		bps = 16;
		break;
	case SHORT:
		bps = 16;
		break;
	case FLOAT:
		bps = 32;
		break;
	default:
		itkExceptionMacro(
			<< "TIFF supports unsigned/signed char, unsigned/signed short, and float");
	}

	uint16_t predictor;

	const char *mode = "w";

	// If the size of the image is greater then 2GB then use big tiff
	const SizeType oneKiloByte = 1024;
	const SizeType oneMegaByte = 1024 * oneKiloByte;
	const SizeType oneGigaByte = 1024 * oneMegaByte;
	const SizeType twoGigaBytes = 2 * oneGigaByte;

	if (this->GetImageSizeInBytes() > twoGigaBytes)
	{
#ifdef TIFF_INT64_T  // detect if libtiff4
		// Adding the "8" option enables the use of big tiff
		mode = "w8";
#else
		itkExceptionMacro(<< "Size of image exceeds the limit of libtiff.");
#endif
	}

	TIFF *tif = TIFFOpen(m_FileName.c_str(), mode);
	if (!tif)
	{
		itkExceptionMacro("Error while trying to open file for writing: "
			<< this->GetFileName()
			<< std::endl
			<< "Reason: "
			<< itksys::SystemTools::GetLastSystemError());
	}

	if (this->GetComponentType() == SHORT
		|| this->GetComponentType() == CHAR)
	{
		TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
	}
	else if (this->GetComponentType() == FLOAT)
	{
		TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
	}

	uint32 w = width;
	uint32 h = height;

	if (m_NumberOfDimensions == 3)
	{
		TIFFCreateDirectory(tif);
	}
	for (page = 0; page < pages; page++)
	{
		TIFFSetDirectory(tif, page);
		TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w);
		TIFFSetField(tif, TIFFTAG_IMAGELENGTH, h);
		TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
		TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, scomponents);
		TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bps); // Fix for stype
		TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
		if (this->GetComponentType() == SHORT
			|| this->GetComponentType() == CHAR)
		{
			TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
		}
		else if (this->GetComponentType() == FLOAT)
		{
			TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
		}
		TIFFSetField(tif, TIFFTAG_SOFTWARE, "InsightToolkit");

		if (scomponents > 3)
		{
			// if number of scalar components is greater than 3, that means we assume
			// there is alpha.
			uint16  extra_samples = scomponents - 3;
			uint16 *sample_info = new uint16[scomponents - 3];
			sample_info[0] = EXTRASAMPLE_ASSOCALPHA;
			int cc;
			for (cc = 1; cc < scomponents - 3; cc++)
			{
				sample_info[cc] = EXTRASAMPLE_UNSPECIFIED;
			}
			TIFFSetField(tif, TIFFTAG_EXTRASAMPLES, extra_samples,
				sample_info);
			delete[] sample_info;
		}

		int compression;

		if (m_UseCompression)
		{
			switch (m_Compression)
			{
			case TIFFImageIO::LZW:
				itkWarningMacro(<< "LZW compression is patented outside US so it is disabled. packbits compression will be used instead");
				ITK_FALLTHROUGH;
			case TIFFImageIO::PackBits:
				compression = COMPRESSION_PACKBITS; break;
			case TIFFImageIO::JPEG:
				compression = COMPRESSION_JPEG; break;
			case TIFFImageIO::Deflate:
				compression = COMPRESSION_DEFLATE; break;
			default:
				compression = COMPRESSION_NONE;
			}
		}
		else
		{
			compression = COMPRESSION_NONE;
		}

		TIFFSetField(tif, TIFFTAG_COMPRESSION, compression); // Fix for compression

		uint16 photometric = (scomponents == 1) ? PHOTOMETRIC_MINISBLACK : PHOTOMETRIC_RGB;

		if (compression == COMPRESSION_JPEG)
		{
			TIFFSetField(tif, TIFFTAG_JPEGQUALITY, m_JPEGQuality);
			TIFFSetField(tif, TIFFTAG_JPEGCOLORMODE, JPEGCOLORMODE_RGB);
		}
		else if (compression == COMPRESSION_DEFLATE)
		{
			predictor = 2;
			TIFFSetField(tif, TIFFTAG_PREDICTOR, predictor);
		}

		TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, photometric); // Fix for scomponents

		// Previously, rowsperstrip was set to a default value so that it would be calculated using
		// the STRIP_SIZE_DEFAULT defined to be 8 kB in tiffiop.h.
		// However, this a very conservative small number, and it leads to very small strips resulting
		// in many io operations, which can be slow when written over networks that require
		// encryption/decryption of each packet (such as sshfs).
		// Conversely, if the value is too high, a lot of extra memory is required to store the strips
		// before they are written out.
		// Experiments writing TIFF images to drives mapped by sshfs showed that a good tradeoff is
		// achieved when the STRIP_SIZE_DEFAULT is increased to 1 MB.
		// This results in an increase in memory usage but no increase in writing time when writing
		// locally and significant writing time improvement when writing over sshfs.
		// For example, writing a 2048x2048 uint16 image with 8 kB per strip leads to 2 rows per strip
		// and takes about 120 seconds writing over sshfs.
		// Using 1 MB per strip leads to 256 rows per strip, which takes only 4 seconds to write over sshfs.
		// Rather than change that value in the third party libtiff library, we instead compute the
		// rowsperstrip here to lead to this same value.
#ifdef TIFF_INT64_T // detect if libtiff4
		uint64_t scanlinesize = TIFFScanlineSize64(tif);
#else
		tsize_t scanlinesize = TIFFScanlineSize(tif);
#endif
		if (scanlinesize == 0)
		{
			itkExceptionMacro("TIFFScanlineSize returned 0");
		}
		rowsperstrip = (uint32_t)(1024 * 1024 / scanlinesize);
		if (rowsperstrip < 1)
		{
			rowsperstrip = 1;
		}

		TIFFSetField(tif,
			TIFFTAG_ROWSPERSTRIP,
			TIFFDefaultStripSize(tif, rowsperstrip));

		if (resolution_x > 0 && resolution_y > 0)
		{
			TIFFSetField(tif, TIFFTAG_XRESOLUTION, resolution_x);
			TIFFSetField(tif, TIFFTAG_YRESOLUTION, resolution_y);
			TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH);
		}

		if (m_NumberOfDimensions == 3)
		{
			// We are writing single page of the multipage file
			TIFFSetField(tif, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
			// Set the page number
			TIFFSetField(tif, TIFFTAG_PAGENUMBER, page, pages);
		}
		int rowLength; // in bytes

		switch (this->GetComponentType())
		{
		case UCHAR:
			rowLength = sizeof(unsigned char);
			break;
		case USHORT:
			rowLength = sizeof(unsigned short);
			break;
		case CHAR:
			rowLength = sizeof(char);
			break;
		case SHORT:
			rowLength = sizeof(short);
			break;
		case FLOAT:
			rowLength = sizeof(float);
			break;
		default:
			itkExceptionMacro(
				<< "TIFF supports unsigned/signed char, unsigned/signed short, and float");
		}

		rowLength *= this->GetNumberOfComponents();
		rowLength *= width;

		int row = 0;
		for (unsigned int idx2 = 0; idx2 < height; idx2++)
		{
			if (TIFFWriteScanline(tif, const_cast<char *>(outPtr), row, 0) < 0)
			{
				itkExceptionMacro(<< "TIFFImageIO: error out of disk space");
			}
			outPtr += rowLength;
			++row;
		}

		if (m_NumberOfDimensions == 3)
		{
			TIFFWriteDirectory(tif);
		}
	}
	TIFFClose(tif);
}


wTIFFImageIOFactory::wTIFFImageIOFactory()
{
	// std::cout << "wTIFFImageIOFactory::wTIFFImageIOFactory()" << std::endl;

	this->RegisterOverride("itkImageIOBase",
		"wTiffIO",
		"TIFF Image IO (Unicode)",
		1,
		itk::CreateObjectFunction< wTIFFImageIO >::New());
}

wTIFFImageIOFactory::~wTIFFImageIOFactory()
{}

const char *
wTIFFImageIOFactory::GetITKSourceVersion(void) const
{
	return ITK_SOURCE_VERSION;
}

const char *
wTIFFImageIOFactory::GetDescription(void) const
{
	return "TIFF ImageIO Factory (Unicode), allows the loading of TIFF images into insight";
}

// Automatically register wTIFFImageIOFactory with ITK on startup,
// making wTIFFImageIO selectable for ImageFileReader.
static class wTiffSelfRegistration
{
public:
	wTiffSelfRegistration()
	{
		// Find and un-register the original factory
		std::string oldname = "TIFFImageIOFactory";
		std::list<itk::ObjectFactoryBase *> factories = itk::ObjectFactoryBase::GetRegisteredFactories();
		for (itk::ObjectFactoryBase *f : factories) {
			if (oldname == f->GetNameOfClass()) {
				// std::cout << "Unregister " << f->GetNameOfClass() << std::endl;
				itk::ObjectFactoryBase::UnRegisterFactory(f);
				break;
			}
		}
		// Now register the unicode version
		wTIFFImageIOFactory::RegisterOneFactory();
	}
} tiffSelfRegistration;

#endif
