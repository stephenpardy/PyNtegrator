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

        int err;
	PyArrayObject *mass_gal,
	              *rad_gal,
		      *pm_mu_delta,
		      *pm_mu_alphacosdelta,
		      *distance_gal,
		      *l,
		      *b;
	double tpast,
	       sigma_x,
	       sigma_v,
	       sigma_vx,
	       sigma_mu;

	// Parse the input tuple
	// Change theser arguments to accept python lists
	if (!PyArg_ParseTuple(args, "OOOOOOOddddd",
			      &mass_gal,
			      &rad_gal,
			      &pm_mu_delta,
			      &pm_mu_alphacosdelta,
			      &distance_gal,
			      &l,
			      &b,
			      &tpast,
			      &sigma_x,
			      &sigma_v,
			      &sigma_vx,
			      &sigma_mu))	// reads in input parameters
		return NULL;

	// Call the external C function to create the streakline model and return the likelihood
	
	err = orbit(mass_gal,
	            rad_gal,
		    pm_mu_delta,
		    pm_mu_alphacosdelta,
		    distance_gal,
		    l,
		    b,
		    tpast,
		    sigma_x,
		    sigma_v,
		    sigma_vx,
		    sigma_mu);

    // Check if error raised
	if(err!=0) {
		PyErr_SetString(PyExc_RuntimeError, "Error occured in the leapfrog integrator.");
		return NULL;
	}

	// Return likelihood value
	return Py_BuildValue("d", 0.0);

}
