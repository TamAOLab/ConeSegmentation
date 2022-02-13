#ifdef _WIN32

#include <sys/stat.h>
#include <locale>
#include <codecvt>
#include <itkVersion.h>
#include <itksys/SystemTools.hxx>
#include <itkByteSwapper.h>
#include <iostream>

#include "wBmpIO.h"

void wBMPImageIO::OpenFileForReading(std::ifstream & inputStream, const std::string & filename,
	bool ascii)
{
	// Make sure that we have a file to
	if (filename.empty())
	{
		itkExceptionMacro(<< "A FileName must be specified.");
	}

	// Close file from any previous image
	if (inputStream.is_open())
	{
		inputStream.close();
	}

	// Open the new file for reading
	itkDebugMacro(<< "Opening file for reading: " << filename);

	std::ios::openmode mode = std::ios::in;
	if (!ascii)
	{
		mode |= std::ios::binary;
	}

	std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
	std::string fn(filename);
	std::wstring wfn = converter.from_bytes(fn);
	inputStream.open(wfn, mode);

	if (!inputStream.is_open() || inputStream.fail())
	{
		itkExceptionMacro(<< "Could not open file: "
			<< filename << " for reading."
			<< std::endl
			<< "Reason: "
			<< itksys::SystemTools::GetLastSystemError());
	}
}

void wBMPImageIO::OpenFileForWriting(std::ofstream & outputStream, const std::string & filename,
	bool truncate, bool ascii)
{
	// Make sure that we have a file to
	if (filename.empty())
	{
		itkExceptionMacro(<< "A FileName must be specified.");
	}

	// Close file from any previous image
	if (outputStream.is_open())
	{
		outputStream.close();
	}

	// Open the new file for writing
	itkDebugMacro(<< "Opening file for writing: " << filename);

	std::ios::openmode mode = std::ios::out;
	if (truncate)
	{
		// typically, ios::out also implies ios::trunc, but being explicit is safer
		mode |= std::ios::trunc;
	}
	else
	{
		mode |= std::ios::in;
		// opening a nonexistent file for reading + writing is not allowed on some platforms
		if (!itksys::SystemTools::FileExists(filename.c_str()))
		{
			itksys::SystemTools::Touch(filename.c_str(), true);
			// don't worry about failure here, errors should be detected later when the file
			// is "actually" opened, unless there is a race condition
		}
	}
	if (!ascii)
	{
		mode |= std::ios::binary;
	}

	std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
	std::string fn(filename);
	std::wstring wfn = converter.from_bytes(fn);
	outputStream.open(wfn, mode);

	if (!outputStream.is_open() || outputStream.fail())
	{
		itkExceptionMacro(<< "Could not open file: "
			<< filename << " for writing."
			<< std::endl
			<< "Reason: "
			<< itksys::SystemTools::GetLastSystemError());
	}
}


wBMPImageIOFactory::wBMPImageIOFactory()
{
	this->RegisterOverride("itkImageIOBase",
		"wBmpIO",
		"BMP Image IO (Unicode)",
		1,
		itk::CreateObjectFunction< wBMPImageIO >::New());
}

wBMPImageIOFactory::~wBMPImageIOFactory()
{}

const char *
wBMPImageIOFactory::GetITKSourceVersion(void) const
{
	return ITK_SOURCE_VERSION;
}

const char *
wBMPImageIOFactory::GetDescription(void) const
{
	return "BMP ImageIO Factory (Unicode), allows the loading of BMP images into insight";
}

// Automatically register wBMPImageIOFactory with ITK on startup,
// making wBMPImageIO selectable for ImageFileReader.
static class wBmpSelfRegistration
{
public:
	wBmpSelfRegistration()
	{
		// Find and un-register the original factory
		std::string oldname = "BMPImageIOFactory";
		std::list<itk::ObjectFactoryBase *> factories = itk::ObjectFactoryBase::GetRegisteredFactories();
		for (itk::ObjectFactoryBase *f : factories) {
			if (oldname == f->GetNameOfClass()) {
				// std::cout << "Unregister " << f->GetNameOfClass() << std::endl;
				itk::ObjectFactoryBase::UnRegisterFactory(f);
				break;
			}
		}
		// Now register the unicode version
		wBMPImageIOFactory::RegisterOneFactory();
	}
} bmpSelfRegistration;

#endif
