#include <Python.h>
#include "python.h"

void dCopy_Array(double* array, int len, PyObject* list)
{
  int i;
  Py_ssize_t pylen=len;
  //check before copy
  if(PyList_Size(list)<pylen){
    printf("Error while copying double C array into Python list: not enough place!\n");
  }
  else{
    for(i=0; i<len; i++){
      PyList_SetItem(list, i, PyFloat_FromDouble(array[i]));
    }
  }
}

void dInit_Array(PyObject* list)
{
  int i;
  for(i=0; i<PyList_Size(list); i++){
    PyList_SetItem(list,i,PyFloat_FromDouble(0));
  }
}

void dPaste_Array(double* array, int len, PyObject* list)
{
  int i;
  Py_ssize_t pylen=len;
  Py_ssize_t pyi;
  //check before copy
  if(PyList_Size(list)<pylen){
    printf("Error while copying double C array into Python list: not enough place!\n");
  }
  else{
    for(i=0; i<len; i++){ //I don't trust pyis
      pyi=i;
      array[i]=PyFloat_AsDouble(PyList_GetItem(list,pyi));
    }
  }
}

void eigvals(double* A, int dim, double* wr, double* wi,
	     double* vecl, double* vecr)
{
  Py_ssize_t Pydim=dim;
  Py_ssize_t Pydimdim=dim*dim;
  PyObject* PyA=PyList_New(Pydimdim);
  PyObject* Pywr=PyList_New(Pydim);
  PyObject* Pywi=PyList_New(Pydim);
  PyObject* Pyvecl=PyList_New(Pydimdim);
  PyObject* Pyvecr=PyList_New(Pydimdim);
  PyObject* main_module;
  PyObject* global_dict;
  PyObject* expression;
  PyObject* pArgs;
  PyObject* pNewArgs;
  PyObject* Keys;
  int i;

  dCopy_Array(A, dim*dim, PyA);
  dInit_Array(Pywr);
  dInit_Array(Pywi);
  dInit_Array(Pyvecl);
  dInit_Array(Pyvecr);
  //Do python stuff
  main_module= PyImport_AddModule("__main__");
  global_dict = PyModule_GetDict(main_module);
  Keys=PyDict_Keys(global_dict);
  expression = PyDict_GetItemString(global_dict, "eigvals");
  if(expression==NULL){
    printf("Could not load function\n");
  }
  pArgs= PyTuple_Pack(5,PyA,Pywr,Pywi,Pyvecl,Pyvecr);
  if(pArgs==NULL){
    printf("This is NULL!\n");
  }
  pNewArgs=PyObject_CallObject(expression,pArgs);
  if(pNewArgs==NULL){
    printf("This is NULL!\n");
  }
  PyA=PyTuple_GetItem(pNewArgs,0);
  Pywr=PyTuple_GetItem(pNewArgs,1);
  Pywi=PyTuple_GetItem(pNewArgs,2);
  Pyvecl=PyTuple_GetItem(pNewArgs,3);
  Pyvecr=PyTuple_GetItem(pNewArgs,4);
  Py_DECREF(pArgs);

  //then copy back
  dPaste_Array(wr, dim, Pywr);
  dPaste_Array(wi, dim, Pywi);
  dPaste_Array(vecl, dim*dim, Pyvecl);
  dPaste_Array(vecr, dim*dim, Pyvecr);
  Py_DECREF(PyA);
  Py_DECREF(Pywr);
  Py_DECREF(Pyvecl);
  Py_DECREF(Pyvecr);
  Py_DECREF(Pywi);
}
