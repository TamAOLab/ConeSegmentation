#ifndef wbmpio_h
#define wbmpio_h
#ifdef _WIN32

#include <itkBMPImageIO.h>
#include <itkObjectFactoryBase.h>

class wBMPImageIO : public itk::BMPImageIO
{
public:
	/** Standard class typedefs. */
	typedef wBMPImageIO                Self;
	typedef itk::BMPImageIO            Superclass;
	typedef itk::SmartPointer< Self >  Pointer;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(BMPImageIO, Superclass);

protected:
	virtual void OpenFileForReading(std::ifstream & inputStream, const std::string & filename,
		bool ascii = false) override;
	virtual void OpenFileForWriting(std::ofstream & outputStream, const std::string & filename,
		bool truncate = true, bool ascii = false) override;
};


class wBMPImageIOFactory : public itk::ObjectFactoryBase
{
public:
	/** Standard class typedefs. */
	typedef wBMPImageIOFactory         Self;
	typedef itk::ObjectFactoryBase          Superclass;
	typedef itk::SmartPointer< Self >       Pointer;
	typedef itk::SmartPointer< const Self > ConstPointer;

	/** Class methods used to interface with the registered factories. */
	virtual const char * GetITKSourceVersion(void) const ITK_OVERRIDE;

	virtual const char * GetDescription(void) const ITK_OVERRIDE;

	/** Method for class instantiation. */
	itkFactorylessNewMacro(Self);
	static wBMPImageIOFactory * FactoryNew() { return new wBMPImageIOFactory; }
	/** Run-time type information (and related methods). */
	itkTypeMacro(wBMPImageIOFactory, ObjectFactoryBase);

	/** Register one factory of this type  */
	static void RegisterOneFactory(void)
	{
		wBMPImageIOFactory::Pointer BMPFactory = wBMPImageIOFactory::New();

		itk::ObjectFactoryBase::RegisterFactoryInternal(BMPFactory);
	}

protected:
	wBMPImageIOFactory();
	~wBMPImageIOFactory() ITK_OVERRIDE;

private:
	ITK_DISALLOW_COPY_AND_ASSIGN(wBMPImageIOFactory);
};

#endif
#endif
