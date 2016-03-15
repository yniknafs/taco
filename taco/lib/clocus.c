//
// TACO
// clocus.c
//
#include <Python.h>

#include "cgtf.h"


static PyObject *
py_gtf_index_loci(PyObject *self, PyObject *args) {
    char *gtf_file;
    char *output_file;
    int ret;

    if (!PyArg_ParseTuple(args, "ss", &gtf_file, &output_file)) {
        return NULL;
    }
    ret = gtf_index_loci(gtf_file, output_file);
    return Py_BuildValue("i", ret);
}

static PyMethodDef CLocusMethods[] = {
    {"gtf_index_loci", py_gtf_index_loci, METH_VARARGS},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initclocus(void)
{
    (void) Py_InitModule("clocus", CLocusMethods);
}
