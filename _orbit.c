#include <Python.h>
#include <numpy/arrayobject.h>
/* Docstrings */
static char module_docstring[] = "This module provides an interface for running a Runge-Kutta streakline integrator using C.";
static char orbit_docstring[] = "Generate streakline model in various potentials and compare them to observational data.";

/* Available functions */
static PyObject *orbit_orbit(PyObject *self, PyObject *args);

/* Module specification */
static PyMethodDef module_methods[] = {
	{"orbit", orbit_orbit, METH_VARARGS, orbit_docstring},
	{NULL, NULL, 0, NULL}
};

/* Initialize the module */
PyMODINIT_FUNC init_orbit(void)
{
	PyObject *m = Py_InitModule3("_orbit", module_methods, module_docstring);
	if (m == NULL)
		return;

	/* Load `numpy` functionality. */
	import_array();
}

static PyObject *orbit_orbit(PyObject *self, PyObject *args)
{

        int err, input_type, int_mode, ngals;
        PyDictObject *parameters;
       	// Parse the input tuple
	// Change theser arguments to accept python lists
	if (!PyArg_ParseTuple(args, "iiiO", &input_type, &int_mode, &ngals, &parameters)) // reads in input parameters
		return NULL;
	 	
        double* output_pos = (double *) malloc(ngals*3*sizeof(double));
	double* output_vel = (double *) malloc(ngals*3*sizeof(double));
        // Call the external C function to create the streakline model and return the likelihood
	err = orbit(int_mode, ngals, parameters, output_pos, output_vel);

    // Check if error raised
	if(err!=0) {
		PyErr_SetString(PyExc_RuntimeError, "Error occured in the leapfrog integrator.");
		return NULL;
	}
        int n, i;
        int array_len = 3*ngals;
        PyObject *lst_pos = PyList_New(array_len);
        PyObject *lst_vel = PyList_New(array_len);
        if (!lst_pos)
            return NULL;
        if (!lst_vel)
            return NULL;

        for (n=0; n<ngals; n++){
            for(i=0; i<3; i++){
                // First get positions.
                PyObject *posnum = PyFloat_FromDouble(output_pos[i+n*3]);
                if (!posnum) {
                    Py_DECREF(lst_pos);
                    return NULL;
                }                  
                PyList_SET_ITEM(lst_pos, i+n*3, posnum);   // reference to num stolen  
                //Now get velocities
                PyObject *velnum = PyFloat_FromDouble(output_vel[i+n*3]);
                if (!velnum) {
                    Py_DECREF(lst_vel);
                    return NULL;
                }
                PyList_SET_ITEM(lst_vel, i+n*3, velnum);   // reference to num stolen  
            }
        }
        PyDictObject *Outputs = PyDict_New();
        PyDict_SetItemString(Outputs, "pos", lst_pos);
        PyDict_SetItemString(Outputs, "vel", lst_vel);
        return Outputs;
	// Return likelihood value
	//return Py_BuildValue("d", 0.0);

}
