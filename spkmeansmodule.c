#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/arrayobject.h>

#include "spkmeans.h"
#include "safelib.h"

static PyObject *py_WAM(PyObject *self, PyObject *args)
{
  npy_intp dimensions[2];
  PyObject *input_datapoints;
  PyArrayObject *res;
  PyArrayObject *datapoints;
  npy_intp i = 0, j = 0;
  matrix_t *datapoints_c = (matrix_t *)safeMalloc(sizeof(matrix_t));
  matrix_t *res_c = (matrix_t *)safeMalloc(sizeof(matrix_t));

  if (!PyArg_ParseTuple(args, "OII:wam_wrapper", &input_datapoints, &datapoints_c->rows, &datapoints_c->cols))
  {
    return NULL;
  }

  datapoints_c->values = (double *)safeMalloc(datapoints_c->rows * datapoints_c->cols * sizeof(double));

  dimensions[0] = datapoints_c->rows;
  dimensions[1] = datapoints_c->rows;

  datapoints = (PyArrayObject *)PyArray_ContiguousFromObject(input_datapoints, NPY_DOUBLE, 2, 2);

  if ((datapoints == NULL))
  {
    return NULL;
  }

  for (i = 0; i < datapoints_c->rows; i++)
  {
    for (j = 0; j < datapoints_c->cols; j++)
    {
      datapoints_c->values[((int)i) * datapoints_c->cols + (int)j] = *(double *)PyArray_GETPTR2(datapoints, i, j);
    }
  }

  WAM(datapoints_c, res_c);

  res = (PyArrayObject *)PyArray_SimpleNew(2, dimensions, NPY_DOUBLE);
  memcpy(PyArray_DATA(res), res_c->values, sizeof(double) * sqr(res_c->rows));

  Py_DECREF(datapoints);
  free(datapoints_c->values);
  free(datapoints_c);
  free(res_c->values);
  free(res_c);

  return (PyObject *)res;
}

static PyObject *py_DDG(PyObject *self, PyObject *args)
{
  npy_intp dimensions[2];
  matrix_t *datapoints_c = (matrix_t *)safeMalloc(sizeof(matrix_t));
  matrix_t *wam_c = (matrix_t *)safeMalloc(sizeof(matrix_t));
  double *ddg_c, *res_c;
  PyObject *input_datapoints;
  PyArrayObject *res;
  PyArrayObject *datapoints;
  npy_intp i = 0, j = 0;

  if (!PyArg_ParseTuple(args, "OII:wam_wrapper", &input_datapoints, &datapoints_c->rows, &datapoints_c->cols))
  {
    return NULL;
  }

  datapoints_c->values = (double *)safeMalloc(datapoints_c->rows * datapoints_c->cols * sizeof(double));
  wam_c->values = (double *)safeCalloc(datapoints_c->rows, datapoints_c->rows * sizeof(double));
  res_c = (double *)safeCalloc(datapoints_c->rows, datapoints_c->rows * sizeof(double));
  ddg_c = (double *)safeCalloc(datapoints_c->rows, sizeof(double));

  dimensions[0] = datapoints_c->rows;
  dimensions[1] = datapoints_c->rows;

  datapoints = (PyArrayObject *)PyArray_ContiguousFromObject(input_datapoints, NPY_DOUBLE, 2, 2);

  if ((datapoints == NULL))
  {
    return NULL;
  }

  for (i = 0; i < datapoints_c->rows; i++)
  {
    for (j = 0; j < datapoints_c->cols; j++)
    {
      datapoints_c->values[((int)i) * datapoints_c->cols + (int)j] = *(double *)PyArray_GETPTR2(datapoints, i, j);
    }
  }

  WAM(datapoints_c, wam_c);
  DDG(wam_c, ddg_c);

  for (i = 0; i < wam_c->rows; i++)
  {
    res_c[i * wam_c->cols + i] = ddg_c[i];
  }

  res = (PyArrayObject *)PyArray_SimpleNew(2, dimensions, NPY_DOUBLE);
  memcpy(PyArray_DATA(res), res_c, sizeof(double) * sqr(wam_c->rows));

  Py_DECREF(datapoints);
  free(datapoints_c->values);
  free(datapoints_c);
  free(wam_c->values);
  free(wam_c);
  free(ddg_c);
  free(res_c);

  return (PyObject *)res;
}

static PyObject *py_lNorm(PyObject *self, PyObject *args)
{
  npy_intp dimensions[2];
  matrix_t *datapoints_c = (matrix_t *)safeMalloc(sizeof(matrix_t));
  matrix_t *wam_c = (matrix_t *)safeMalloc(sizeof(matrix_t));
  matrix_t *res_c = (matrix_t *)safeMalloc(sizeof(matrix_t));
  double *ddg_c;
  PyObject *input_datapoints;
  PyArrayObject *res;
  PyArrayObject *datapoints;
  npy_intp i = 0, j = 0;

  if (!PyArg_ParseTuple(args, "OII:wam_wrapper", &input_datapoints, &datapoints_c->rows, &datapoints_c->cols))
  {
    return NULL;
  }

  datapoints_c->values = (double *)safeMalloc(datapoints_c->rows * datapoints_c->cols * sizeof(double));
  wam_c->values = (double *)safeCalloc(datapoints_c->rows, datapoints_c->rows * sizeof(double));
  res_c->values = (double *)safeCalloc(datapoints_c->rows, datapoints_c->rows * sizeof(double));
  ddg_c = (double *)safeCalloc(datapoints_c->rows, sizeof(double));

  dimensions[0] = datapoints_c->rows;
  dimensions[1] = datapoints_c->rows;

  datapoints = (PyArrayObject *)PyArray_ContiguousFromObject(input_datapoints, NPY_DOUBLE, 2, 2);

  if ((datapoints == NULL))
  {
    return NULL;
  }

  for (i = 0; i < datapoints_c->rows; i++)
  {
    for (j = 0; j < datapoints_c->cols; j++)
    {
      datapoints_c->values[((int)i) * datapoints_c->cols + (int)j] = *(double *)PyArray_GETPTR2(datapoints, i, j);
    }
  }

  WAM(datapoints_c, wam_c);
  DDG(wam_c, ddg_c);
  DHalf(ddg_c, datapoints_c->rows);
  laplacian(wam_c, ddg_c, res_c);

  res = (PyArrayObject *)PyArray_SimpleNew(2, dimensions, NPY_DOUBLE);
  memcpy(PyArray_DATA(res), res_c->values, sizeof(double) * sqr(res_c->rows));

  Py_DECREF(datapoints);
  free(datapoints_c->values);
  free(datapoints_c);
  free(wam_c->values);
  free(wam_c);
  free(res_c->values);
  free(res_c);
  free(ddg_c);

  return (PyObject *)res;
}

static PyObject *py_fit(PyObject *self, PyObject *args)
{
  uint32_t datasetSize;
  uint32_t dim;
  uint32_t clusterCount;
  npy_intp dimensions[2];
  double *centroids_c;
  double *datapoints_c;
  PyObject *input_centroids, *input_datapoints;
  PyArrayObject *res;
  PyArrayObject *centroids, *datapoints;
  npy_intp i = 0, j = 0;

  if (!PyArg_ParseTuple(args, "OOIII:fit_wrapper", &input_centroids, &input_datapoints, &datasetSize, &dim, &clusterCount))
  {
    return NULL;
  }

  centroids_c = (double *)safeCalloc(clusterCount, dim * sizeof(double));
  datapoints_c = (double *)safeCalloc(datasetSize, dim * sizeof(double));

  dimensions[0] = clusterCount;
  dimensions[1] = dim;

  centroids = (PyArrayObject *)PyArray_ContiguousFromObject(input_centroids, NPY_DOUBLE, 2, 2);
  datapoints = (PyArrayObject *)PyArray_ContiguousFromObject(input_datapoints, NPY_DOUBLE, 2, 2);

  if ((centroids == NULL) || (datapoints == NULL))
  {
    return NULL;
  }

  for (i = 0; i < datasetSize; i++)
  {
    for (j = 0; j < dim; j++)
    {
      datapoints_c[((int)i) * dim + (int)j] = *(double *)PyArray_GETPTR2(datapoints, i, j);
    }
  }

  for (i = 0; i < clusterCount; i++)
  {
    for (j = 0; j < dim; j++)
    {
      centroids_c[((int)i) * dim + (int)j] = *(double *)PyArray_GETPTR2(centroids, i, j);
    }
  }

  kmeansFit(centroids_c, datapoints_c, datasetSize, dim, clusterCount);

  res = (PyArrayObject *)PyArray_SimpleNew(2, dimensions, NPY_DOUBLE);
  memcpy(PyArray_DATA(res), centroids_c, sizeof(double) * clusterCount * dim);

  Py_DECREF(centroids);
  Py_DECREF(datapoints);
  free(datapoints_c);
  free(centroids_c);

  return (PyObject *)res;
}

static PyObject *py_jacobi(PyObject *self, PyObject *args)
{
  npy_intp dimensions[2];
  npy_intp eigenCount[1];
  double *V_c, *eigenArray_c;
  matrix_t *matrix_c = (matrix_t *)safeMalloc(sizeof(matrix_t));
  PyObject *input_matrix;
  PyArrayObject *V;
  PyArrayObject *matrix;
  PyArrayObject *eigenArray;
  npy_intp i = 0, j = 0;

  if (!PyArg_ParseTuple(args, "OI:wam_wrapper", &input_matrix, &matrix_c->rows))
  {
    return NULL;
  }

  matrix_c->cols = matrix_c->rows;

  matrix_c->values = (double *)safeCalloc(matrix_c->rows, matrix_c->cols * sizeof(double));
  V_c = (double *)safeCalloc(matrix_c->rows, matrix_c->cols * sizeof(double));
  eigenArray_c = (double *)safeCalloc(matrix_c->rows, sizeof(double));

  dimensions[0] = matrix_c->rows;
  dimensions[1] = matrix_c->rows;

  matrix = (PyArrayObject *)PyArray_ContiguousFromObject(input_matrix, NPY_DOUBLE, 2, 2);

  if ((matrix == NULL))
  {
    return NULL;
  }

  for (i = 0; i < matrix_c->rows; i++)
  {
    for (j = 0; j < matrix_c->cols; j++)
    {
      matrix_c->values[((int)i) * matrix_c->cols + (int)j] = *(double *)PyArray_GETPTR2(matrix, i, j);
    }
  }

  Jacobi(matrix_c, V_c, eigenArray_c);

  eigenCount[0] = (npy_intp)matrix_c->rows;
  eigenArray = (PyArrayObject *)PyArray_SimpleNew(1, eigenCount, NPY_DOUBLE);
  V = (PyArrayObject *)PyArray_SimpleNew(2, dimensions, NPY_DOUBLE);
  memcpy(PyArray_DATA(eigenArray), eigenArray_c, sizeof(double) * matrix_c->rows);
  memcpy(PyArray_DATA(V), V_c, sizeof(double) * sqr(matrix_c->rows));

  Py_DECREF(matrix);
  free(eigenArray_c);
  free(matrix_c->values);
  free(matrix_c);
  free(V_c);

  return Py_BuildValue("OO", eigenArray, V);
}

static PyObject *py_prepate_data(PyObject *self, PyObject *args)
{
  uint32_t k;
  npy_intp dimensions[2];
  matrix_t *datapoints_c = (matrix_t *)safeMalloc(sizeof(matrix_t));
  matrix_t *res_c = (matrix_t *)safeMalloc(sizeof(matrix_t));
  PyObject *input_datapoints;
  PyArrayObject *res;
  PyArrayObject *datapoints;
  npy_intp i = 0, j = 0;

  if (!PyArg_ParseTuple(args, "OIII:wam_wrapper", &input_datapoints, &datapoints_c->rows, &datapoints_c->cols, &k))
  {
    return NULL;
  }

  datapoints_c->values = (double *)safeCalloc(datapoints_c->rows, datapoints_c->cols * sizeof(double));

  datapoints = (PyArrayObject *)PyArray_ContiguousFromObject(input_datapoints, NPY_DOUBLE, 2, 2);

  if ((datapoints == NULL))
  {
    return NULL;
  }

  for (i = 0; i < datapoints_c->rows; i++)
  {
    for (j = 0; j < datapoints_c->cols; j++)
    {
      datapoints_c->values[((int)i) * datapoints_c->cols + (int)j] = *(double *)PyArray_GETPTR2(datapoints, i, j);
    }
  }

  res_c = prepareData(datapoints_c, k);

  dimensions[0] = res_c->rows;
  dimensions[1] = res_c->cols;

  res = (PyArrayObject *)PyArray_SimpleNew(2, dimensions, NPY_DOUBLE);
  memcpy(PyArray_DATA(res), res_c->values, sizeof(double) * res_c->rows * res_c->cols);

  Py_DECREF(datapoints);

  free(datapoints_c->values);
  free(datapoints_c);
  free(res_c->values);
  free(res_c);

  return (PyObject *)res;
}

static PyMethodDef capiMethods[] = {
    {"fit",
     (PyCFunction)py_fit,
     METH_VARARGS,
     PyDoc_STR("A function used to fit, performing the kmeans algorithm")},
    {"WAM",
     (PyCFunction)py_WAM,
     METH_VARARGS,
     PyDoc_STR("A function used to calculate the WAM from a set of observations")},
    {"DDG",
     (PyCFunction)py_DDG,
     METH_VARARGS,
     PyDoc_STR("A function used to calculate the DDG from a set of observations")},
    {"lNorm",
     (PyCFunction)py_lNorm,
     METH_VARARGS,
     PyDoc_STR("A function used to calculate the lNorm of the WAM from a set of observations")},
    {"jacobi",
     (PyCFunction)py_jacobi,
     METH_VARARGS,
     PyDoc_STR("A function used to calculate the eigen values and vectors of a given matrix")},
    {"prepare_data",
     (PyCFunction)py_prepate_data,
     METH_VARARGS,
     PyDoc_STR("A function used to prepare data from set of observations")},
    {NULL, NULL, 0, NULL}};

/* This initiates the module using the above definitions. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeans", /* name of module */
    NULL,
    -1,
    capiMethods};

PyMODINIT_FUNC PyInit_spkmeans(void)
{
  PyObject *m;

  m = PyModule_Create(&moduledef);
  if (!m)
  {
    return NULL;
  }
  import_array();
  return m;
}
