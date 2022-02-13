#ifndef wpngio_h
#define wpngio_h
#ifdef _WIN32

#include <itkPNGImageIO.h>
#include <itkObjectFactoryBase.h>

class wPNGImageIO : public itk::PNGImageIO
{
public:
	typedef wPNGImageIO				Self;
	typedef PNGImageIO				Superclass;
	typedef itk::SmartPointer< Self >	Pointer;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(wPNGImageIO, PNGImageIO);

	virtual bool CanReadFile(const char *file) override;
	virtual void Read(void *buffer) override;
	virtual void ReadImageInformation() override;
	virtual void Write(const void *buffer) override;
protected:
	void WriteSlice(const std::string & fileName, const void *buffer);
};

class wPNGImageIOFactory : public itk::ObjectFactoryBase
{
public:
	/** Standard class typedefs. */
	typedef wPNGImageIOFactory         Self;
	typedef itk::ObjectFactoryBase          Superclass;
	typedef itk::SmartPointer< Self >       Pointer;
	typedef itk::SmartPointer< const Self > ConstPointer;

	/** Class methods used to interface with the registered factories. */
	virtual const char * GetITKSourceVersion(void) const ITK_OVERRIDE;

	virtual const char * GetDescription(void) const ITK_OVERRIDE;

	/** Method for class instantiation. */
	itkFactorylessNewMacro(Self);
	static wPNGImageIOFactory * FactoryNew() { return new wPNGImageIOFactory; }
	/** Run-time type information (and related methods). */
	itkTypeMacro(wPNGImageIOFactory, ObjectFactoryBase);

	/** Register one factory of this type  */
	static void RegisterOneFactory(void)
	{
		wPNGImageIOFactory::Pointer PNGFactory = wPNGImageIOFactory::New();

		itk::ObjectFactoryBase::RegisterFactoryInternal(PNGFactory);
	}

protected:
	wPNGImageIOFactory();
	~wPNGImageIOFactory() ITK_OVERRIDE;

private:
	ITK_DISALLOW_COPY_AND_ASSIGN(wPNGImageIOFactory);
};

#endif
#endif
