%{
#include <array>
#include <type_traits>

template<typename T, std::size_t N>
static bool CheckTuple(PyObject * input, const bool writeErrorMessage) {
  if (!PySequence_Check(input)) {
    if (writeErrorMessage)
      PyErr_SetString(PyExc_TypeError, "Expecting a sequence");
    return false;
  }
  if (PyObject_Length(input) != N) {
    if (writeErrorMessage)
      PyErr_SetString(PyExc_ValueError, "Sequence size mismatch");
    return false;
  }
  for (unsigned int i = 0; i < N; ++i) {
    PyObject *o = PySequence_GetItem(input, i);

    if (std::is_same<T, double>::value) {
      if (!PyFloat_Check(o)) {
        Py_XDECREF(o);
        if (writeErrorMessage)
          PyErr_SetString(PyExc_ValueError, "Expecting a sequence of floats");
        return false;
      }
      //array->at(i) = PyFloat_AsDouble(o);
      Py_DECREF(o);
    } else if (std::is_same<T, int>::value) {
      if (!PyInt_Check(o)) {
        Py_XDECREF(o);
        if (writeErrorMessage)
          PyErr_SetString(PyExc_ValueError, "Expecting a sequence of integers");
        return false;
      }
      //array->at(i) = PyInt_AsLong(o);
      Py_DECREF(o);
    } else {
      Py_XDECREF(o);
      if (writeErrorMessage)
        PyErr_SetString(PyExc_ValueError, "unknown type");
      return false;
    }
  }
  return true;
}

template<typename T, std::size_t N>
static bool TupleToStdArray(PyObject * input, std::array<T, N> * const array) {
  if (!CheckTuple<T, N>(input, true))
    return false;

  for (unsigned int i = 0; i < N; ++i) {
    PyObject *o = PySequence_GetItem(input, i);

    if (std::is_same<T, double>::value) {
      array->at(i) = PyFloat_AsDouble(o);
      Py_DECREF(o);
    } else if (std::is_same<T, int>::value) {
      array->at(i) = PyInt_AsLong(o);
      Py_DECREF(o);
    }
  }
  return true;
}

template<typename T, std::size_t N>
static PyObject * StdArrayToTuple(const std::array<T, N>& array) {
  PyObject *o = PyTuple_New(N);
  if (std::is_same<T, int>::value) {
    for (unsigned int i = 0; i < N; ++i) {
      PyTuple_SetItem(o, i, PyInt_FromLong(array[i]));
    }
  } else if (std::is_same<T, double>::value) {
    for (unsigned int i = 0; i < N; ++i) {
      PyTuple_SetItem(o, i, PyFloat_FromDouble(array[i]));
    }
  } else {
    PyErr_SetString(PyExc_ValueError, "unknown type");
    return NULL;
  }
  return o;
}

%}

%define %std_array_typemap(type, size)
%typemap(in) std::array<type, size> {
  if (!TupleToStdArray<type, size>($input, &$1))
    return NULL;
}

%typemap(in) std::array<type, size>& (std::array<type, size> temp) {
  if (!TupleToStdArray<type, size>($input, &temp))
    return NULL;
  $1 = &temp;
}

%typemap(varin) std::array<type, size> {
  if (!TupleToStdArray<type, size>($input, &$1))
    return NULL;
}

%typemap(varin) std::array<type, size>& (std::array<type, size> temp) {
  if (!TupleToStdArray<type, size>($input, &temp))
    return NULL;
  $1 = &temp;
}

%typemap(typecheck) std::array<type, size> {
  $1 = CheckTuple<type, size>($input, false) ? 1 : 0;
}

%typemap(typecheck) std::array<type, size>& (std::array<type, size> temp) {
  $1 = CheckTuple<type, size>($input, false) ? 1 : 0;
}

%typemap(out) std::array<type, size> {
  $result = StdArrayToTuple<type, size>($1);
}

%typemap(out) std::array<type, size>& {
  $result = StdArrayToTuple<type, size>(*$1);
}

%typemap(varout) std::array<type, size> {
  $result = StdArrayToTuple<type, size>($1);
}

%typemap(varout) std::array<type, size>& {
  $result = StdArrayToTuple<type, size>(*$1);
}
%enddef

namespace std {
template<typename T, std::size_t N>
struct array {};

%extend array {
  inline size_t __len__() const { return $self->size(); }

  inline const T& __getitem__(size_t i) const throw(std::out_of_range) {
    if (i >= N)
      throw std::out_of_range("out of bounds access");
    return $self->data()[i];
  }

  inline void __setitem__(size_t i, const T& v) throw(std::out_of_range) {
    if (i >= N)
      throw std::out_of_range("out of bounds access");
    $self->data()[i] = v;
  }
}

}
