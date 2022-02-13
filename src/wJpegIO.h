#ifndef wjpegio_h
#define wjpegio_h
#ifdef _WIN32

#include <itkJPEGImageIO.h>
#include <itkObjectFactoryBase.h>

class wJPEGImageIO :public itk::JPEGImageIO
{
public:
	typedef wJPEGImageIO          Self;
	typedef JPEGImageIO          Superclass;
	typedef itk::SmartPointer< Self > Pointer;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro(wJPEGImageIO, JPEGImageIO);

	virtual bool CanReadFile(const char *file) override;
	virtual void Read(void *buffer) override;
	virtual void ReadImageInformation() override;
	virtual void Write(const void *buffer) override;
protected:
	void WriteSlice(std::string & fileName, const void *buffer);
};

class wJPEGImageIOFactory : public itk::ObjectFactoryBase
{
public:
	/** Standard class typedefs. */
	typedef wJPEGImageIOFactory         Self;
	typedef itk::ObjectFactoryBase          Superclass;
	typedef itk::SmartPointer< Self >       Pointer;
	typedef itk::SmartPointer< const Self > ConstPointer;

	/** Class methods used to interface with the registered factories. */
	virtual const char * GetITKSourceVersion(void) const ITK_OVERRIDE;

	virtual const char * GetDescription(void) const ITK_OVERRIDE;

	/** Method for class instantiation. */
	itkFactorylessNewMacro(Self);
	static wJPEGImageIOFactory * FactoryNew() { return new wJPEGImageIOFactory; }
	/** Run-time type information (and related methods). */
	itkTypeMacro(wJPEGImageIOFactory, ObjectFactoryBase);

	/** Register one factory of this type  */
	static void RegisterOneFactory(void)
	{
		wJPEGImageIOFactory::Pointer JPEGFactory = wJPEGImageIOFactory::New();

		itk::ObjectFactoryBase::RegisterFactoryInternal(JPEGFactory);
	}

protected:
	wJPEGImageIOFactory();
	~wJPEGImageIOFactory() ITK_OVERRIDE;

private:
	ITK_DISALLOW_COPY_AND_ASSIGN(wJPEGImageIOFactory);
};


#endif
#endif
