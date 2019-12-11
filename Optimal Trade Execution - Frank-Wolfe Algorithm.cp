#include <gurobi_c++.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

double gradient(double* vector, double alpha, double pi, int T) {
	double retcode = 0, sum = 0;
	int i;
	double* v;
	double* x = NULL, * impact = NULL, * temp = NULL, * temp0 = NULL;

	/*store the value of vector to x*/
	x = (double*)calloc(2 * T, sizeof(double));

	for (int i = 0; i < T; i++) {
		x[i] = vector[i];
	}
	//std::cout << "x address" << x << "\n x value" << *x;


	/*calcuate the impact fuction of vector*/
	for (i = 0; i < T; i++) {
		vector[i] = 1 - alpha * pow(vector[i], pi);
	}
	/*store the value of impact fuction to impact*/
	impact = (double*)calloc(2 * T, sizeof(double));

	for (int i = 0; i < T; i++) {
		impact[i] = vector[i];
	}


	/*calcuate the cumprod of the impact fuction of vector*/
	v = vector;

	for (i = 1; i < T; i++) {
		v[i] = vector[i] * v[i - 1];
	}

	/*store the cumprod to temp0*/
	temp0 = (double*)calloc(2 * T, sizeof(double));

	for (int i = 0; i < T; i++) {
		temp0[i] = v[i];
	}

	/*calculate the intermediate term */
	for (i = 0; i < T; i++) {
		v[i] = v[i] * x[i];
	}

	/*store the intermediate term to temp*/
	temp = (double*)calloc(2 * T, sizeof(double));

	for (int i = 0; i < T; i++) {
		temp[i] = v[i];
	}

	/*calcuate term1 */
	for (i = 0; i < T; i++) sum += v[i];
	v[0] = sum;
	for (i = 1; i < T; i++) {
		v[i] = v[i - 1] - temp[i - 1];
	}
	for (i = 0; i < T; i++) {
		v[i] = v[i] / impact[i];
	}

	/*calcuate deri_impact */
	for (i = 0; i < T; i++) {
		x[i] = -alpha * pi * pow(x[i], pi - 1);
	}
	/*calcuate part2 */
	for (i = 0; i < T; i++) {
		v[i] = v[i] * x[i];
	}
	/*calcuate gradient to v*/
	for (i = 0; i < T; i++) {
		v[i] = temp0[i] + v[i];
	}


	/*for (int i = 0; i < T; i++) {
		printf("gradient = %.20f", v[i]);
		printf("\n");
	}*/


	return retcode;
}

double my_function(double* vector, double alpha, double pi, int T) {
	double fun_val = 0;
	int i;
	double* v;
	double* x = NULL;

	/*store the value of vector to x*/
	x = (double*)calloc(2 * T, sizeof(double));

	for (int i = 0; i < T; i++) {
		x[i] = vector[i];
	}

	double* vector_copy = NULL;
	vector_copy = (double*)calloc(2 * T, sizeof(double));

	for (int i = 0; i < T; i++) {
		vector_copy[i] = vector[i];
	}

	/*calcuate the impact fuction of vector*/
	for (i = 0; i < T; i++) {
		vector_copy[i] = 1 - alpha * pow(vector_copy[i], pi);
	}

	/*calcuate the cumprod of the impact fuction of vector*/

	for (i = 1; i < T; i++) {
		vector_copy[i] = vector_copy[i] * vector_copy[i - 1];
	}

	/*calcuate the dop product of cumprod and x*/
	for (i = 0; i < T; i++) {
		fun_val += vector_copy[i] * x[i];
	}
	return fun_val;
}

double compute_slacks(double* vector, int N, int T) {
	double slacks = N;
	int i;
	for (i = 0; i < T; i++) {
		slacks -= vector[i];
	}

	return slacks;
}


double line_search(double* vector, double* y, double alpha, double pi, int T, double step_size) {

	double max_t = 0, f = 0, t = 0;
	int i, j;
	double* x_new = NULL;
	x_new = (double*)calloc(2 * T, sizeof(double));



	for (double i = 0; i <= 1; i = i + step_size) {

		/*x is the new vector after moving*/
		for (int j = 0; j < T; j++) {
			x_new[j] = vector[j] + i * y[j];
		}

		f = my_function(x_new, alpha, pi, T);

		if (f > max_t) {
			max_t = f;
			t = i;
		}
	}

	return t;
}


int direction(double* grad_vec, double* y, int T, int N, double* vector)
{
	int retcode = 0;
	GRBenv* env = NULL;
	GRBmodel* model = NULL;
	int n, j, Nq;
	double* obj = NULL;
	double* lb = NULL;
	double* ub = NULL;
	double* x;
	int* cind, * qrow, * qcol;
	double rhs;
	char sense;
	double* cval, * qval;
	int numnonz;

	n = T;
	retcode = GRBloadenv(&env, "first.log");
	if (retcode) goto BACK;

	/* Create initial model */
	retcode = GRBnewmodel(env, &model, "first", n,
		NULL, NULL, NULL, NULL, NULL);
	if (retcode) goto BACK;

	/** next we create the 4 columns **/
	obj = (double*)calloc(n, sizeof(double));
	ub = (double*)calloc(n, sizeof(double));
	lb = (double*)calloc(n, sizeof(double));
	x = (double*)calloc(n, sizeof(double));

	for (j = 0; j < n; j++) {
		obj[j] = -grad_vec[j];               //min -gradient
	}
	//obj[0] = 1; obj[1] = 2; obj[2] = 3; obj[3] = 4; obj[4] = 5;

	for (j = 0; j < n; j++) {
		ub[j] = N - vector[j];
	}
	for (j = 0; j < n; j++) {
		lb[j] = -vector[j];
	}

	/* initialize variables */
	for (j = 0; j < n; j++) {
		retcode = GRBsetstrattrelement(model, "VarName", j, "NewConstr");
		if (retcode) goto BACK;

		retcode = GRBsetdblattrelement(model, "Obj", j, obj[j]);
		if (retcode) goto BACK;

		retcode = GRBsetdblattrelement(model, "LB", j, lb[j]);
		if (retcode) goto BACK;

		retcode = GRBsetdblattrelement(model, "UB", j, ub[j]);
		if (retcode) goto BACK;
	}

	/** now we will add one constraint at a time **/
	/** we need to have a couple of auxiliary arrays **/

	cind = (int*)calloc(n, sizeof(int));
	cval = (double*)calloc(n, sizeof(double));


	/** constraint is sum all x = 0 **/
	for (j = 0; j < n; j++) {
		cind[j] = j;
	}
	for (j = 0; j < n; j++) {
		cval[j] = 1;
	}

	numnonz = n;
	rhs = 0;         // which equal the slacks
	sense = GRB_EQUAL;

	retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, "first_constraint");
	if (retcode) goto BACK;

	retcode = GRBupdatemodel(model);
	if (retcode) goto BACK;


	/** optional: write the problem **/

	//retcode = GRBwrite(model, "myfirst.lp");
	//if (retcode) goto BACK;


	retcode = GRBoptimize(model);
	if (retcode) goto BACK;


	/** get solution **/


	retcode = GRBgetdblattrarray(model,
		GRB_DBL_ATTR_X, 0, n,
		x);
	if (retcode) goto BACK;

	/** now let's see the values **/

	//for (j = 0; j < n; j++) {
	//	printf("variables = %g\n", x[j]);
	//}


	GRBfreeenv(env);

	for (int i = 0; i < T; i++) {
		y[i] = x[i];
	}


BACK:
	printf("\nexiting with retcode %d\n", retcode);
	return retcode;
}


int main() {
	int code = 0, slacks = 0;
	double alpha = 0.0001, pi = 0.5, c;
	double* vector = NULL;
	int i, N = 10000, T = 20;
	double* v;
	double step_size = 1E-4;
	int time = 1;

	/*initialize feasible solution x.(called vector)*/
	vector = (double*)calloc(2 * T, sizeof(double));


	for (int i = 0; i < T; i++) {
		c = (double)N / T;
		vector[i] = c;
	}

	//Get the original function value
	double old_f = my_function(vector, alpha, pi, T);  //vector1 NOT changed now;
	printf("old funciton:%f", old_f);
	std::cout << "This is the " << time << " iteration";
	time++;


	//Get the gradient for the current position vector
	double* grad_vec = NULL;
	grad_vec = (double*)calloc(2 * T, sizeof(double));

	for (int i = 0; i < T; i++) {
		grad_vec[i] = vector[i];
	}

	gradient(grad_vec, alpha, pi, T);   // vector becomes the gradient

	//for (int i = 0; i < T; i++) {
	//	printf("gradient = %.20f", grad_vec[i]);
	//	printf("\n");
	//}

	slacks = compute_slacks(grad_vec, N, T);    // which always equals 0. So we are not using the computed slack in optimization. (just set to 0)


	/*compute step direction using gurobi, changing y */
	double* y = NULL;
	y = (double*)calloc(T, sizeof(double));
	direction(grad_vec, y, T, N, vector);

	// y is the step direction now;
	for (int i = 0; i < T; i++) {
		std::cout << "y " << i << " = " << y[i] << "\n";
	}

	//line search to get the optimal t step
	double t = line_search(vector, y, alpha, pi, T, step_size);
	printf("t = %f", t);
	std::cout << "\n";


	//update the original position to the new position
	for (int i = 0; i < T; i++) {
		vector[i] += t * y[i];
		std::cout << "The Amount to be sold on day " << i << ":   " << vector[i] << "\n";
	}

	//get new function value after the
	double new_f = my_function(vector, alpha, pi, T);
	printf("new_function = %f", new_f);
	std::cout << "\n";

	double diff = new_f - old_f;
	std::cout << "The Value Improvement is: " << diff << "\n";


	//Continue iteration until reach certain condition
	while (diff > 1E-10) {

		std::cout << "This is the " << time << " iteration" << "\n";
		time++;

		//Get the gradient for the current position vector
		double* grad_vec = NULL;
		grad_vec = (double*)calloc(2 * T, sizeof(double));

		for (int i = 0; i < T; i++) {
			grad_vec[i] = vector[i];
		}

		gradient(grad_vec, alpha, pi, T);   // vector becomes the gradient

		//for (int i = 0; i < T; i++) {
		//	printf("gradient = %.20f", grad_vec[i]);
		//	printf("\n");
		//}

		slacks = compute_slacks(grad_vec, N, T);    // which always equals 0. So we are not using the computed slack in optimization. (just set to 0)

		/*compute step direction using gurobi, changing y */
		double* y = NULL;
		y = (double*)calloc(T, sizeof(double));
		direction(grad_vec, y, T, N, vector);

		// y is the step direction now;
		for (int i = 0; i < T; i++) {
			std::cout << "y " << i << " = " << y[i] << "\n";
		}

		//line search to get the optimal t step
		double t = line_search(vector, y, alpha, pi, T, step_size);
		printf("t = %f", t);
		std::cout << "\n";


		//update the original position to the new position
		for (int i = 0; i < T; i++) {
			vector[i] += t * y[i];
			std::cout << "The Amount to be sold on day " << i << ":   " << vector[i] << "\n";
		}

		//avoid deep copy by updating diff first
		diff = my_function(vector, alpha, pi, T) - new_f;

		new_f = my_function(vector, alpha, pi, T);

		printf("new_function value = %f", new_f);
		std::cout << "\n";
		std::cout << "The Value Improvement is: " << diff << "\n";

	}

	return code;
}


