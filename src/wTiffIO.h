#ifndef wtiffio_h
#define wtiffio_h
#ifdef _WIN32

#include <itkTIFFImageIO.h>
#include <itkObjectFactoryBase.h>

class wTIFFImageIO : public itk::TIFFImageIO
{
public:
	typedef wTIFFImageIO          Self;
	typedef TIFFImageIO          Superclass;
	typedef itk::SmartPointer< Self > Pointer;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(wTIFFImageIO, TIFFImageIO);

	virtual bool CanReadFile(const char *file) override;
	virtual void Write(const void *buffer) override;

protected:
	void InternalWrite(const void *buffer);
};

class wTIFFImageIOFactory :public itk::ObjectFactoryBase
{
public:
	/** Standard class typedefs. */
	typedef wTIFFImageIOFactory         Self;
	typedef itk::ObjectFactoryBase          Superclass;
	typedef itk::SmartPointer< Self >       Pointer;
	typedef itk::SmartPointer< const Self > ConstPointer;

	/** Class methods used to interface with the registered factories. */
	virtual const char * GetITKSourceVersion(void) const ITK_OVERRIDE;

	virtual const char * GetDescription(void) const ITK_OVERRIDE;

	/** Method for class instantiation. */
	itkFactorylessNewMacro(Self);
	static wTIFFImageIOFactory * FactoryNew() { return new wTIFFImageIOFactory; }
	/** Run-time type information (and related methods). */
	itkTypeMacro(wTIFFImageIOFactory, ObjectFactoryBase);

	/** Register one factory of this type  */
	static void RegisterOneFactory(void)
	{
		wTIFFImageIOFactory::Pointer TIFFFactory = wTIFFImageIOFactory::New();

		itk::ObjectFactoryBase::RegisterFactoryInternal(TIFFFactory);
	}

protected:
	wTIFFImageIOFactory();
	~wTIFFImageIOFactory() ITK_OVERRIDE;

private:
	ITK_DISALLOW_COPY_AND_ASSIGN(wTIFFImageIOFactory);
};


#endif
#endif