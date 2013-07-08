#ifndef FITS_READER_H
#define FITS_READER_H

#include <vtkImageReader2.h>
#include <vtkPVConfig.h>

class FitsReader : public vtkImageReader2
{
public:
    static FitsReader* New();
    vtkTypeMacro(FitsReader,vtkImageReader2);
    virtual void PrintSelf(ostream& os, vtkIndent indent);

    // Description:
    // Can the file be loaded with this plugin
    virtual int CanReadFile(const char* fname);

    // Description:
    // Get the file extensions for this format.
    // Returns a string with a space separated list of extensions in
    // the format .extension
    virtual const char* GetFileExtensions()
    {
        return ".fits .hdf5";
    }

    // Description:
    // Return a descriptive name for the file format that might be useful in a GUI.
    virtual const char* GetDescriptiveName()
    {
        return "FITS";
    }

protected:
    FitsReader();
    virtual ~FitsReader();

    virtual void ExecuteInformation();
#if PARAVIEW_VERSION_MINOR <= 20
    virtual void ExecuteData(vtkDataObject* out);
#else
    virtual void ExecuteDataWithInformation(vtkDataObject* data, vtkInformation* outInfo);
#endif
private:
    FitsReader(const FitsReader&);  // Not implemented.
    void operator=(const FitsReader&);  // Not implemented.
};
#endif // FITS_READER_H
