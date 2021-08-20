#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/ndarraytypes.h>
#include <numpy/arrayobject.h>

#include "spkmeans.h"

static PyObject *py_WAM(PyObject *self, PyObject *args)
{
    uint32_t count;
    uint32_t dim;
    npy_intp dimensions[2];
    double *datapoints_c, *res_c;
    PyObject *input_datapoints;
    PyArrayObject *res;
    PyArrayObject *datapoints;
    npy_intp i = 0, j = 0;

    if (!PyArg_ParseTuple(args, "OII:wam_wrapper", &input_datapoints, &count, &dim))
    {
        return NULL;
    }

    datapoints_c = (double *) calloc(count, dim * sizeof(double));
    res_c = (double *) calloc(count, count * sizeof(double));

    assert(datapoints_c != NULL && res_c != NULL);

    dimensions[0] = count;
    dimensions[1] = count;

    datapoints = (PyArrayObject *)PyArray_ContiguousFromObject(input_datapoints, NPY_DOUBLE, 2, 2);

    if ((datapoints == NULL))
    {
        return NULL;
    }

    for (i = 0; i < count; i++)
    {
        for (j = 0; j < dim; j++)
        {
            datapoints_c[((int)i) * dim + (int)j] = *(double *)PyArray_GETPTR2(datapoints, i, j);
        }
    }

    WAM(datapoints_c, count, dim, res_c);

    res = (PyArrayObject *)PyArray_SimpleNew(2, dimensions, NPY_DOUBLE);
    memcpy(PyArray_DATA(res), res_c, sizeof(double) * sqr(count));

    Py_DECREF(datapoints);
    free(datapoints_c);
    free(res_c);

    return (PyObject *)res;
}

static PyObject *py_DDG(PyObject *self, PyObject *args)
{
    uint32_t count;
    uint32_t dim;
    npy_intp dimensions[2];
    double *datapoints_c, *wam_c, *ddg_c, *res_c;
    PyObject *input_datapoints;
    PyArrayObject *res;
    PyArrayObject *datapoints;
    npy_intp i = 0, j = 0;

    if (!PyArg_ParseTuple(args, "OII:wam_wrapper", &input_datapoints, &count, &dim))
    {
        return NULL;
    }

    datapoints_c = (double *) calloc(count, dim * sizeof(double));
    wam_c = (double *)calloc(count, count * sizeof(double));
    ddg_c = (double *) calloc(count, sizeof(double));

    assert(datapoints_c != NULL && wam_c != NULL && ddg_c != NULL);

    dimensions[0] = count;
    dimensions[1] = count;

    datapoints = (PyArrayObject *)PyArray_ContiguousFromObject(input_datapoints, NPY_DOUBLE, 2, 2);

    if ((datapoints == NULL))
    {
        return NULL;
    }

    for (i = 0; i < count; i++)
    {
        for (j = 0; j < dim; j++)
        {
            datapoints_c[((int)i) * dim + (int)j] = *(double *)PyArray_GETPTR2(datapoints, i, j);
        }
    }

    WAM(datapoints, count, dim, wam_c);
    DDG(wam_c, count, ddg_c);
    res_c = prettyDDG(ddg_c, count);

    res = (PyArrayObject *)PyArray_SimpleNew(2, dimensions, NPY_DOUBLE);
    memcpy(PyArray_DATA(res), res_c, sizeof(double) * sqr(count));

    Py_DECREF(datapoints);
    free(datapoints_c);
    free(wam_c);
    free(ddg_c);
    free(res_c);

    return (PyObject *)res;
}

static PyObject *py_lNorm(PyObject *self, PyObject *args)
{
    uint32_t count;
    uint32_t dim;
    npy_intp dimensions[2];
    double *datapoints_c, *wam_c, *ddg_c, *res_c;
    PyObject *input_datapoints;
    PyArrayObject *res;
    PyArrayObject *datapoints;
    npy_intp i = 0, j = 0;

    if (!PyArg_ParseTuple(args, "OII:wam_wrapper", &input_datapoints, &count, &dim))
    {
        return NULL;
    }

    datapoints_c = (double *) calloc(count, dim * sizeof(double));
    wam_c = (double *)calloc(count, count * sizeof(double));
    ddg_c = (double *) calloc(count, sizeof(double));

    assert(datapoints_c != NULL && wam_c != NULL && ddg_c != NULL);

    dimensions[0] = count;
    dimensions[1] = count;

    datapoints = (PyArrayObject *)PyArray_ContiguousFromObject(input_datapoints, NPY_DOUBLE, 2, 2);

    if ((datapoints == NULL))
    {
        return NULL;
    }

    for (i = 0; i < count; i++)
    {
        for (j = 0; j < dim; j++)
        {
            datapoints_c[((int)i) * dim + (int)j] = *(double *)PyArray_GETPTR2(datapoints, i, j);
        }
    }

    WAM(datapoints, count, dim, wam_c);
    DDG(wam_c, count, ddg_c);
    DHalf(ddg_c, count);
    laplacian(wam_c, ddg_c, count, res_c);

    res = (PyArrayObject *)PyArray_SimpleNew(2, dimensions, NPY_DOUBLE);
    memcpy(PyArray_DATA(res), res_c, sizeof(double) * sqr(count));

    Py_DECREF(datapoints);
    free(datapoints_c);
    free(wam_c);
    free(ddg_c);
    free(res_c);

    return (PyObject *)res;
}

static PyObject *py_fit(PyObject *self, PyObject *args)
{
    uint32_t MAX_ITER;
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

    if (!PyArg_ParseTuple(args, "OOIIII:fit_wrapper", &input_centroids, &input_datapoints, &MAX_ITER, &datasetSize, &dim, &clusterCount))
    {
        return NULL;
    }

    centroids_c = calloc(clusterCount, dim * sizeof(double));
    datapoints_c = calloc(datasetSize, dim * sizeof(double));

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

    kmeansFit(centroids_c, datapoints_c, datasetSize, dim, clusterCount, MAX_ITER);

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
    uint32_t count;
    npy_intp dimensions[2];
    double *matrix_c, *V_c, *eigenArray_c;
    PyObject *input_matrix;
    PyArrayObject *V;
    PyArrayObject *matrix;
    PyArrayObject *eigenArray;
    npy_intp i = 0, j = 0;

    if (!PyArg_ParseTuple(args, "OI:wam_wrapper", &input_matrix, &count))
    {
        return NULL;
    }

    matrix_c = (double *) calloc(count, count * sizeof(double));
    V_c = (double *) calloc(count, count * sizeof(double));
    eigenArray_c = (double *) calloc(count, sizeof(double));

    assert(matrix_c != NULL && V_c != NULL && eigenArray_c != NULL);

    dimensions[0] = count;
    dimensions[1] = count;

    matrix = (PyArrayObject *)PyArray_ContiguousFromObject(input_matrix, NPY_DOUBLE, 2, 2);

    if ((matrix == NULL))
    {
        return NULL;
    }

    for (i = 0; i < count; i++)
    {
        for (j = 0; j < count; j++)
        {
            matrix_c[((int)i) * count + (int)j] = *(double *)PyArray_GETPTR2(matrix, i, j);
        }
    }

    Jacobi(matrix_c, count, V_c, eigenArray);

    eigenArray = (PyArrayObject *)PyArray_SimpleNew(1, (npy_int *){count}, NPY_DOUBLE);
    V = (PyArrayObject *)PyArray_SimpleNew(2, dimensions, NPY_DOUBLE);
    memcpy(PyArray_DATA(eigenArray), eigenArray_c, sizeof(double) * count);
    memcpy(PyArray_DATA(V), V_c, sizeof(double) * sqr(count));

    Py_DECREF(matrix);
    free(eigenArray_c);
    free(matrix_c);
    free(V_c);

    return Py_BuildValue("OO",eigenArray,V);
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
    {"Jacobi",
     (PyCFunction)py_jacobi,
     METH_VARARGS,
     PyDoc_STR("A function used to calculate the eigen values and vectors of a given matrix")},
    {NULL, NULL, 0, NULL}};

/* This initiates the module using the above definitions. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeans", /* name of module */
    NULL,
    -1,
    capiMethods};

PyMODINIT_FUNC PyInit_mykmeanspp(void)
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
