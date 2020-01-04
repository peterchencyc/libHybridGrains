// PythonObject.cpp
//
// Breannan Smith
// Last updated: 09/22/2015

#include "PythonObject.h"

#ifdef USE_PYTHON

PythonObject::PythonObject(PyObject *object) : m_object(object) {}

PythonObject::PythonObject(PythonObject &&other) : m_object(other.m_object) {
  other.m_object = nullptr;
}

PythonObject::~PythonObject() { Py_XDECREF(m_object); }

PythonObject &PythonObject::operator=(PythonObject &&other) {
  if (this != &other) {
    // We are discarding this object's PyObject so we must decrease the
    // PyObject's reference count
    Py_XDECREF(m_object);
    m_object = other.m_object;
    other.m_object = nullptr;
  }
  return *this;
}

PythonObject::operator PyObject *() const { return m_object; }

#endif
