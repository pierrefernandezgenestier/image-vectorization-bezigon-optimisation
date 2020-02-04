#include "Entities.h"	

using namespace cv;
using namespace std;
using namespace Eigen;

Bezier::Bezier() {
};

Bezier::Bezier(Point2d p0, Point2d p1, Point2d p2, Point2d p3) {
	control_points[0] = p0;	control_points[1] = p1;
	control_points[2] = p2;	control_points[3] = p3;
};

array<Bezier, 2> Bezier::subdivide(double t0) {
	Point2d p0 = control_points[0];
	Point2d p1 = barycenter(t0, p0, control_points[1]);
	Point2d paux = barycenter(t0, control_points[1], control_points[2]);
	Point2d p6 = control_points[3];
	Point2d p5 = barycenter(t0, control_points[2], p6);
	Point2d p2 = barycenter(t0, p1, paux);
	Point2d p4 = barycenter(t0, paux, p5);
	Point2d p3 = barycenter(t0, p2, p4);
	return { Bezier(p0,p1,p2,p3), Bezier(p3,p4,p5,p6) };
}

Bezier Bezier::subdivide(double t0, double t1) {
	return subdivide(t0)[1].subdivide(t1)[1];
};


Point2d Bezier::cubic_interpolation(double t) {
	double s = 1 - t;
	return  s * s * s * control_points[0] + 3 * s * s * t * control_points[1] + 3 * s * t * t * control_points[2] + t * t * t * control_points[3];
}

array<Point2d, 2> Bezier::get_bounding_box() {
	vector<Point2d> pts;
	Point2d c3 = -control_points[0] + 3 * control_points[1] - 3 * control_points[2] + control_points[3];
	Point2d c2 = 3 * control_points[0] - 6 * control_points[1] + 3 * control_points[2];
	Point2d c1 = -3 * control_points[0] + 3 * control_points[1];
	vector<double> roots_dx = get_roots(3 * c3.x, 2 * c2.x, c1.x);
	vector<double> roots_dy = get_roots(3 * c3.y, 2 * c2.y, c1.y);
	for (double root : roots_dx) {
		if (0 <= root && root <= 1) pts.push_back(cubic_interpolation(root));
	}
	for (double root : roots_dy) {
		if (0 <= root && root <= 1) pts.push_back(cubic_interpolation(root));
	}
	pts.push_back(control_points[0]);
	pts.push_back(control_points[3]);
	Point2d corner_min(INFINITY, INFINITY);
	Point2d corner_max(0., 0.);
	for (Point2d pt : pts) {
		if (pt.x < corner_min.x) corner_min.x = pt.x;
		if (pt.y < corner_min.y) corner_min.y = pt.y;
		if (pt.x > corner_max.x) corner_max.x = pt.x;
		if (pt.y > corner_max.y) corner_max.y = pt.y;
	}
	return { corner_min, corner_max };
}

Bezigon::Bezigon() {
}

Bezigon::Bezigon(MatrixXd _Bx, MatrixXd _By, Vec3b _C) {
	Bx = _Bx;
	By = _By;
	lo = get_length();
	C = _C;
};

Bezigon::Bezigon(vector<Point2d> vp, Vec3b _C, int method) {
	if (method == 0) {
		int nb_bezier = vp.size() / 3;
		cout << nb_bezier << endl;
		Bx = MatrixXd(nb_bezier, 3);
		By = MatrixXd(nb_bezier, 3);
		C = _C;
		for (int jj = 0; jj < nb_bezier; jj++) {
			Point2d pt1 = vp[3 * jj + 1];
			Point2d pt2 = vp[3 * jj + 2];
			Bx.row(jj) << vp[3 * jj].x, pt1.x, pt2.x;
			By.row(jj) << vp[3 * jj].y, pt1.y, pt2.y;
		}
		lo = get_length();
	}
	else if (method == 1) {
		int nb_bezier = vp.size();
		Bx = MatrixXd(nb_bezier, 3);
		By = MatrixXd(nb_bezier, 3);
		C = _C;
		for (int jj = 0; jj < nb_bezier - 1; jj++) {
			Point2d pt1 = barycenter(0.33, vp[jj], vp[jj + 1]);
			Point2d pt2 = barycenter(0.66, vp[jj], vp[jj + 1]);
			Bx.row(jj) << vp[jj].x, pt1.x, pt2.x;
			By.row(jj) << vp[jj].y, pt1.y, pt2.y;
		}
		Point2d pt1 = barycenter(0.33, vp[nb_bezier - 1], vp[0]);
		Point2d pt2 = barycenter(0.66, vp[nb_bezier - 1], vp[0]);
		Bx.row(nb_bezier - 1) << vp[nb_bezier - 1].x, pt1.x, pt2.x;
		By.row(nb_bezier - 1) << vp[nb_bezier - 1].y, pt1.y, pt2.y;
		lo = get_length();
	}
	else if (method == 2) {
		int nb_bezier = vp.size();
		Bx = MatrixXd(nb_bezier, 3);
		By = MatrixXd(nb_bezier, 3);
		C = _C;
		Point2f tangent_next, tangent_prev;
		for (int jj = 1; jj < nb_bezier - 1; jj++) {
			tangent_next = vp[jj] + 0.33 * norm(vp[jj + 1] - vp[jj]) * (vp[jj + 1] - vp[jj - 1]) / norm(vp[jj + 1] - vp[jj - 1]);
			tangent_prev = vp[jj] + 0.33 * norm(vp[jj - 1] - vp[jj]) * (vp[jj - 1] - vp[jj + 1]) / norm(vp[jj - 1] - vp[jj + 1]);
			Bx.row(jj - 1)[2] = tangent_prev.x;
			Bx.row(jj)[0] = double(vp[jj].x);
			Bx.row(jj)[1] = tangent_next.x;
			By.row(jj - 1)[2] = tangent_prev.y;
			By.row(jj)[0] = double(vp[jj].y);
			By.row(jj)[1] = tangent_next.y;
		}
		tangent_next = vp[0] + 0.33 * norm(vp[1] - vp[0]) * (vp[1] - vp[nb_bezier - 1]) / norm(vp[1] - vp[nb_bezier - 1]);
		tangent_prev = vp[0] + 0.33 * norm(vp[nb_bezier - 1] - vp[0]) * (vp[nb_bezier - 1] - vp[1]) / norm(vp[1] - vp[nb_bezier - 1]);
		Bx.row(nb_bezier - 1)[2] = tangent_prev.x;
		Bx.row(0)[0] = double(vp[0].x);
		Bx.row(0)[1] = tangent_next.x;
		By.row(nb_bezier - 1)[2] = tangent_prev.y;
		By.row(0)[0] = double(vp[0].y);
		By.row(0)[1] = tangent_next.y;

		tangent_next = vp[nb_bezier - 1] + 0.33 * norm(vp[nb_bezier - 1] - vp[0]) * (vp[0] - vp[nb_bezier - 2]) / norm(vp[0] - vp[nb_bezier - 2]);
		tangent_prev = vp[nb_bezier - 1] + 0.33 * norm(vp[nb_bezier - 1] - vp[nb_bezier - 2]) * (vp[nb_bezier - 2] - vp[0]) / norm(vp[0] - vp[nb_bezier - 2]);
		Bx.row(nb_bezier - 2)[2] = tangent_prev.x;
		Bx.row(nb_bezier - 1)[0] = double(vp[nb_bezier - 1].x);
		Bx.row(nb_bezier - 1)[1] = tangent_next.x;
		By.row(nb_bezier - 2)[2] = tangent_prev.y;
		By.row(nb_bezier - 1)[0] = double(vp[nb_bezier - 1].y);
		By.row(nb_bezier - 1)[1] = tangent_next.y;
		lo = get_length();
	}
}

void Bezigon::set_point_x(int j, int i, double coord_x) {
	Bx.row(j)[i] = coord_x;
};

void Bezigon::set_point_y(int j, int i, double coord_y) {
	By.row(j)[i] = coord_y;
};

Bezier Bezigon::get_bezier(int j) {
	j = (j > 0) ? j % Bx.rows() : (-j) % (-Bx.rows());
	Point2d p0 = Point2d(Bx.row(j)[0], By.row(j)[0]);
	Point2d p1 = Point2d(Bx.row(j)[1], By.row(j)[1]);
	Point2d p2 = Point2d(Bx.row(j)[2], By.row(j)[2]);
	Point2d p3 = (j != Bx.rows() - 1) ?
		Point2d(Bx.row(j + 1)[0], By.row(j + 1)[0]) : Point2d(Bx.row(0)[0], By.row(0)[0]);
	return Bezier(p0, p1, p2, p3);
};

Point2d Bezigon::get_pt(int j, int i = 0) {
	j = (j > 0) ? j % Bx.rows() : (-j) % (-Bx.rows());
	Point2d p0 = Point2d(Bx.row(j)[i], By.row(j)[i]);
	return Point2d(Bx.row(j)[i], By.row(j)[i]);

};

array< Point2d, 2 > Bezigon::get_tangent(int j) {
	Point2d p_prev = get_pt(j - 1, 2);
	Point2d p_1 = get_pt(j, 0);
	Point2d p_2 = get_pt(j, 1);
	return { p_1 - p_prev , p_2 - p_1 };

}

double Bezigon::get_length(int j, double t0, double t1) {
	Bezier bez = get_bezier(j);
	double length = 0.0;
	Point2d previous_point = bez.cubic_interpolation(t0);
	Point2d actual_point;
	for (double t = 0.2; t <= 1.0; t += 0.2) {
		actual_point = bez.cubic_interpolation(t0 + t * (t1 - t0));
		length += distance(actual_point, previous_point);
		previous_point = actual_point;
	}
	return length;
};

double Bezigon::get_length_gl(int j, double t0, double t1) {
	Bezier bez = get_bezier(j);
	double length = 0.;
	Point2d c3 = -bez.control_points[0] + 3 * bez.control_points[1] - 3 * bez.control_points[2] + bez.control_points[3];
	Point2d c2 = 3 * bez.control_points[0] - 6 * bez.control_points[1] + 3 * bez.control_points[2];
	Point2d c1 = -3 * bez.control_points[0] + 3 * bez.control_points[1];
	array<double, 3> ws = { 0.888888, 0.5555555555555556 , 0.5555555555555556 };
	array<double, 3> xs = { 0., -0.7745966692414834 , 0.7745966692414834 };
	double t;
	for (int ii = 0; ii < 3; ii++) {
		t = (t1 - t0) / 2. * xs[ii] + (t1 + t0) / 2.;
		length += (t1 - t0) / 2. * ws[ii] * sqrt(pow(3 * c3.x * t * t + 2 * c2.x * t + c1.x, 2) + pow(3 * c3.y * t * t + 2 * c2.y * t + c1.y, 2));
	}
	return length;
};

double Bezigon::get_length() {
	double length = 0.0;
	for (int jj = 0; jj < Bx.rows(); jj++) {
		//length += get_length(jj);
		length += get_length_gl(jj);
	}
	return length;
};

void Bezigon::update(array<double, 10> vals_inp, int j) {
	if (j < Bx.rows() - 1) {
		set_point_x(j, 1, vals_inp[0]);
		set_point_x(j, 2, vals_inp[1]);
		set_point_x(j + 1, 0, vals_inp[2]);
		set_point_x(j + 1, 1, vals_inp[3]);
		set_point_x(j + 1, 2, vals_inp[4]);

		set_point_y(j, 1, vals_inp[5]);
		set_point_y(j, 2, vals_inp[6]);
		set_point_y(j + 1, 0, vals_inp[7]);
		set_point_y(j + 1, 1, vals_inp[8]);
		set_point_y(j + 1, 2, vals_inp[9]);
	}
	else {
		set_point_x(Bx.rows() - 1, 1, vals_inp[0]);
		set_point_x(Bx.rows() - 1, 2, vals_inp[1]);
		set_point_x(0, 0, vals_inp[2]);
		set_point_x(0, 1, vals_inp[3]);
		set_point_x(0, 2, vals_inp[4]);
		set_point_y(Bx.rows() - 1, 1, vals_inp[5]);
		set_point_y(Bx.rows() - 1, 2, vals_inp[6]);
		set_point_y(0, 0, vals_inp[7]);
		set_point_y(0, 1, vals_inp[8]);
		set_point_y(0, 2, vals_inp[9]);
	}
}

//void Bezigon::rescale(double scale) {
//	Bx = scale * Bx;
//	By = scale * By;
//	lo = scale * lo;
//	// for(int j=0; j <Bx.rows(); j++){
//	// 	for(int i=0;i<3;i++){
//	// 		cv::Point2d point = get_pt(j, i);
//	// 		set_point_x( j,  i, scale*point.x);
//	// 		set_point_y( j,  i, scale*point.y);
//	// 	}
//	// }
//}

array<double, 10> Bezigon::input_propagation(int j) {
	Point2d p1 = get_pt(j, 1);
	Point2d p2 = get_pt(j, 2);
	Point2d p3 = (j != Bx.rows() - 1) ? get_pt(j + 1, 0) : get_pt(0, 0);
	Point2d p4 = (j != Bx.rows() - 1) ? get_pt(j + 1, 1) : get_pt(0, 1);
	Point2d p5 = (j != Bx.rows() - 1) ? get_pt(j + 1, 2) : get_pt(0, 2);
	array<double, 10> x = { double(p1.x),double(p2.x),double(p3.x),double(p4.x),double(p5.x), double(p1.y),double(p2.y),double(p3.y),double(p4.y),double(p5.y) };
	return x;
}

void Bezigon::plot_curve(Image<Vec3b> I, string nom) {
	Bezier bezier_j;
	Image<Vec3b> I_copy;
	I.copyTo(I_copy);
	for (int jj = 0; jj < Bx.rows(); jj++) {
		bezier_j = get_bezier(jj);
		for (double t = 0.0; t <= 1.0; t += 0.05) {
			Point2d m1 = bezier_j.cubic_interpolation(t);
			circle(I_copy, m1, 0, Scalar(0, 255, 0), 2);
			namedWindow(nom, WINDOW_NORMAL);
			resizeWindow(nom, 200, 200);
			imshow(nom, I_copy);
			waitKey(1);
		}
	}
}



VectorizationData::VectorizationData(Bezigon _B, Image<Vec3b> _I) {
	B = _B;
	I = _I;
}



int intersect(Bezier bez1, Bezier bez2, vector<array<double, 2>>* vector_ts, array<double, 2> ts, int n_rec) {
	double thres = 0.1;
	array<Point2d, 2> bbox1 = bez1.get_bounding_box();
	array<Point2d, 2> bbox2 = bez2.get_bounding_box();
	if (overlap(bbox1, bbox2) < 1e-8) return 0;
	else {
		if ((bbox1[1].x - bbox1[0].x) < thres && (bbox1[1].y - bbox1[0].y) < thres) {
			vector_ts->push_back(ts);
			return 1;
		}
		else {
			array<Bezier, 2> bbez1 = bez1.subdivide(0.5);
			array<Bezier, 2> bbez2 = bez2.subdivide(0.5);
			n_rec++;
			return intersect(bbez1[0], bbez2[0], vector_ts, { ts[0],ts[1] }, n_rec)
				+ intersect(bbez1[0], bbez2[1], vector_ts, { ts[0], ts[1] + 1 / pow(2,n_rec) }, n_rec)
				+ intersect(bbez1[1], bbez2[0], vector_ts, { ts[0] + 1 / pow(2,n_rec), ts[1] }, n_rec)
				+ intersect(bbez1[1], bbez2[1], vector_ts, { ts[0] + 1 / pow(2,n_rec), ts[1] + 1 / pow(2,n_rec) }, n_rec);
		}
	}
}

vector<array<double, 2>> intersect(Bezier bez1, Bezier bez2) {
	vector<array<double, 2>> vector_ts;
	int n_int = intersect(bez1, bez2, &vector_ts, { 0.,0. }, 0);
	if (n_int == 0) return vector_ts;

	double eps = 0.05;
	vector<array<double, 2>> final_vector_ts = { vector_ts[0] };
	for (int ii = 1; ii < n_int; ii++) {
		if (abs(vector_ts[ii][0] - vector_ts[ii - 1][0]) > eps&& abs(vector_ts[ii][1] - vector_ts[ii - 1][1]) > eps) final_vector_ts.push_back(vector_ts[ii]);
	}
	return final_vector_ts;
}

vector<array<double, 2>> self_intersect(Bezier bez) {
	vector<array<double, 2>> vector_ts;
	array<Bezier, 2> subdivided = bez.subdivide(0.5);
	int n_int = intersect(subdivided[0], subdivided[1], &vector_ts, { 0.,0.5 }, 1);
	if (n_int == 0) return vector_ts;
	double eps = 0.05;
	vector<array<double, 2>> final_vector_ts = { vector_ts[0] };
	for (int ii = 1; ii < n_int; ii++) {
		if (abs(vector_ts[ii][0] - vector_ts[ii - 1][0]) > eps&& abs(vector_ts[ii][1] - vector_ts[ii - 1][1]) > eps) final_vector_ts.push_back(vector_ts[ii]);
	}
	return final_vector_ts;
}


bool to_svg(vector<VectorizationData> vector_vd, string name) {
	std::ofstream ofs(name);
	if (!ofs.good())
		return false;
	ofs << "<?xml " << "version=\"" << "1.0" << "\" encoding=\"" << "UTF-8" << "\" standalone=\"" << "no" << "\"?>\n "
		<< "<svg " << "xmlns=\"" << "http://www.w3.org/2000/svg" << "\" version=\"" << "1.0" << "\" height=\"" << vector_vd[0].I.height() << "\" width=\"" << vector_vd[0].I.width() << "\" >\n ";
	for (int ii = 0; ii < vector_vd.size(); ii++) {
		VectorizationData vd = vector_vd[ii];
		ofs << "<path id=\"" << ii << "\" \n"
			<< "d=\"M " << vd.B.Bx.row(0)[0] << " " << vd.B.By.row(0)[0] << "\n";
		for (int jj = 0; jj < vd.B.Bx.rows(); jj++) {
			auto cp = vd.B.get_bezier(jj).control_points;
			ofs << "C " << cp[1].x << " " << cp[1].y << " "
				<< cp[2].x << " " << cp[2].y << " "
				<< cp[3].x << " " << cp[3].y << "\n";
		}
		ofs << "\"\n"
			<< "fill=\"rgb(" << (int)(vd.B.C[2]) << "," << (int)vd.B.C[1] << "," << (int)vd.B.C[0] << ")\"\n" << "/>\n";
	}
	ofs << "</svg>";
	ofs.close();
	return true;
}