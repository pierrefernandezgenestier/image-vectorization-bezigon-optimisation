#include "Propagation.h"
#include <ctime>

void propagate(VectorizationData* vd, int nb_iterations, double alpha, double eps, string name) {

	//const double max_update = 0.1;
	/*cout << ">>>Enter Propragation" << std::endl;
	clock_t begin = clock();*/

	////WRITE
	//fstream out;
	//out.open("../test_propagation/" + name + ".txt", fstream::out);
	//out << "Hyperparameters : \n" << "Epsilon : " << eps << "\nAlpha : " << alpha << "\nTaille Image " << vd->I.width() << "x" << vd->I.height();
	//out << "\nLambda_data :" << energy->lambda_data << "\nLambda_handles: " << energy->lambda_handles << "\nLamda_angles :" << energy->lambda_angles << "\nLamda_legnth :" << energy->lambda_length << endl;
	//out << "\nEnergies intialization:";
	//out << "\nE_data :" << energy->energy_data(*vd) << "\nE_handes :" << energy->energy_bezier_handles(vd->B) << "\nE_angles :" << energy->energy_angles(vd->B) << "\nE_length :" << energy->energy_curve_length(vd->B) << endl;;


	double grad, grad_data, grad_prior;
	double engy_plus_data, engy_plus_prior, engy_minus_data, engy_minus_prior;
	array<double, 2> engy_plus, engy_minus;
	array<double, 10> vals_inp, vals_copy_plus, vals_copy_minus, gradients, vals_out;
	//NB OF PROPAGATIONS
	for (int iter = 0; iter < nb_iterations; iter++) {
		//out << "\nPROPAGATION : " << iter << endl;;
		//out<<"\nAlpha_updated :"<<alpha<<endl;;

		//ITERATIONS OVER ALL BEZIERS
		for (int jj = 0; jj < vd->B.Bx.rows(); jj++) {

			//PARAMETERS
			vals_inp = vd->B.input_propagation(jj);
			vals_copy_plus = vals_inp;
			vals_copy_minus = vals_inp;
			//double rescale_max = 1;

			//double engy_prev = energy_partial(*vd, jj);
			//out << "\nBezier Nb: " << jj << "\nPartial_Energy : " << engy_prece << "\nSHIFT : " << endl;
			cout << endl;
			for (int var_index = 0; var_index < 10; var_index++) {

				vals_copy_plus[var_index] = vals_inp[var_index] + eps;
				vals_copy_minus[var_index] = vals_inp[var_index] - eps;

				engy_plus = energy_to_minimize(*vd, jj, vals_copy_plus);
				engy_minus = energy_to_minimize(*vd, jj, vals_copy_minus);

				engy_plus_data = engy_plus[0];
				engy_plus_prior = engy_plus[1];
				engy_minus_data = engy_minus[0];
				engy_minus_prior = engy_minus[1];

				grad_data = (engy_plus_data - engy_minus_data) / (2 * eps);
				grad_prior = (engy_plus_prior - engy_minus_prior) / (2 * eps);

				cout << "grad_prior : " << grad_prior << " , grad_data : " << grad_data << endl;

				grad = grad_data + grad_prior;
				//out << setprecision(6);
				//out << "data= " << grad_data * alpha / double(scale_img) << " prior= " << grad_prior * alpha / double(scale_img);
				//out << " data/prior " << double(grad_data / grad_prior) << endl;

				//if (norm(alpha * grad) / scale_img > max_update) {
				//	if ((norm(alpha * grad) / scale_img) / max_update > rescale_max)
				//		rescale_max = (norm(alpha * grad) / scale_img) / max_update;
				//	//cout<<"update "<<norm(alpha*grad)/scale_img<<endl;
				//	//cout<<"rescale_max "<<rescale_max<<endl;
				//}
				while (norm(alpha * grad) > min(vd->I.width(), vd->I.height()) / 10) {
					grad = grad / 5;
				}
				// cout << "energy_plus: "<<engy_plus_data<<" "<<engy_plus_prior<<endl;	
				// cout << "energy_moins: "<<engy_moins_data<<" "<<engy_moins_prior<<endl;	
				//cout << endl;	

				// cout << endl;
				vals_copy_plus[var_index] = vals_inp[var_index];
				vals_copy_minus[var_index] = vals_inp[var_index];
				gradients[var_index] = grad;
				//cout << "vals_out_before rescale " << vals_inp[var_index] - alpha * grad << endl;
			};
			//cout << endl;
			for (int var_index = 0; var_index < 10; var_index++) {
				//vals_out[var_index] = vals_inp[var_index] - alpha * gradients[var_index] / rescale_max;
				vals_out[var_index] = vals_inp[var_index] - alpha * gradients[var_index];
			}
			vd->B.update(vals_out, jj);
		}
		vd->B.plot_curve(vd->I, " Iter :" + to_string(iter) + " " + name);
		bool svg = to_svg({ *vd }, name + to_string(iter) + ".svg");
	}
	//clock_t end = clock();
	//double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;


	//rms(*vd, " Grande Taille " + name, 5);

	//out << "Propagation time : " << elapsed_secs / 60.0 << " minutes" << endl;
	//out.close();
}