#include "window.h"
#include "func.h"
#include "thread.h"





bool Window::msr_ready() {
	for (int i = 0; i < data.p; i++) {
		if (args[i].ready == false) return false;
	}
	return true;
}

void Window::msr_wait() {
	if (msr_ready()) {
		printf("\n\n%s : Task = 8 R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f It = %d E = %le K = %d Nx = %d Ny = %d P = %d\n", 
			program_name, args[0].r1, args[0].r2, args[0].r3, args[0].r4, args[0].t1, args[0].t2, args[0].its, data.eps, func_id, data.nx, data.ny, data.p);
		//printf("max{|F_min|, |F_max|} = %le\n", f_max);
		/* switch (mode) {
			case 0:		printf("max{|F_min|, |F_max|} = %le\n", F_max); break;
			case 1:		printf("max{|F_min|, |F_max|} = %le\n", PF_max); break;
			case 2:		printf("max{|F_min|, |F_max|} = %le\n", RES_max); break;
		} */
		update();
	}
	else {
		QTimer::singleShot(20, this, &Window::msr_wait);
	}
}

void Window::draw_txt(QPainter *painter) {
	char stroka[100];
	painter->drawText(10, 20, f_name);
	sprintf(stroka, "a  = %.2e", data.a);
	painter->drawText(10, 60, stroka);
	sprintf(stroka, "b  = %.2e", data.b);
	painter->drawText(10, 80, stroka);
	sprintf(stroka, "c  = %.2e", data.c);
	painter->drawText(10, 100, stroka);
	sprintf(stroka, "d  = %.2e", data.d);
	painter->drawText(10, 120, stroka);
	sprintf(stroka, "nx = %d", data.nx);
	painter->drawText(10, 140, stroka);
	sprintf(stroka, "ny = %d", data.ny);
	painter->drawText(10, 160, stroka);
	sprintf(stroka, "mx = %d", data._mx);
	painter->drawText(10, 180, stroka);
	sprintf(stroka, "my = %d", data._my);
	painter->drawText(10, 200, stroka);
	sprintf(stroka, "p  = %d", data.parameter);
	painter->drawText(10, 220, stroka);
}





int Window::parse_command_line(int argc, char *argv[]) {
	if (argc != 13 || sscanf(argv[1], "%lf", &data.a) != 1 || sscanf(argv[2], "%lf", &data.b) != 1 || sscanf(argv[3], "%lf", &data.c) != 1 || sscanf(argv[4], "%lf", &data.d) != 1 || sscanf(argv[5], "%d", &data.nx) != 1 || sscanf(argv[6], "%d", &data.ny) != 1|| sscanf(argv[7], "%d", &data._mx) != 1 || sscanf(argv[8], "%d", &data._my) != 1 || sscanf(argv[9], "%d", &data.func_id) != 1 || sscanf(argv[10], "%lf", &data.eps) != 1 || sscanf(argv[11], "%d", &data.maxit) != 1 || sscanf(argv[12], "%d", &data.p) != 1) return -1;

	if (data._mx <= DEF__MX)	data.mx = data._mx;
	else data.mx = DEF__MX;
	if (data._my <= DEF__MY)	data.my = data._my;
	else data.my = DEF__MY;

	//data._mx = data.mx;	data._my = data.my;

	//data.read_data(argv);
	//printf("epsilon = %lf\n", data.eps);
	if (data.eps < 1e-14) {
		printf("Слишком маленький эпсилон, увеличить его быстро!\n");
		data.eps = DEF_EPS;
	}
	data.repeat_new_data();

	args = new Args[data.p];
	tid = new pthread_t[data.p];
	program_name = argv[0];
	func_id = data.func_id - 1;
	data.func_id = func_id;
	data.parameter = 0;

	//printf("epsilon = %le\n", data.eps);
	//printf("func_id = %d\n", data.func_id);

	// это копируется еще кое-где
	for (int k = 0; k < data.p; k++) {
		args[k].copy_data(data);

		args[k].k = k;
		args[k].mutex = &mutex;
		args[k].cond = &cond;
	}
	
	//one_deletion_of_results();//

	change_func();
    
	return 0;
}

Window::Window(QWidget *parent): QWidget(parent) {
	widget = parent;
}

Window::~Window() {
	if (args) delete[] args;
	if (tid) delete[] tid;
    one_deletion_of_results();
}

QSize Window::minimumSizeHint() const {
	return QSize(100, 100);
}

QSize Window::sizeHint() const {
	return QSize(1000, 1000);
}



void Window::set_f() {
    double (*f)(double, double) = data.f;
    switch(func_id) {
    case 0:
        f_name = "k = 0; f(x, y) = 1";
        f = function_0;
        break;
    case 1:
        f_name = "k = 1; f(x, y) = x";
        f = function_1;
        break;
    case 2:
        f_name = "k = 2; f(x, y) = y";
        f = function_2;
        break;
    case 3:
        f_name = "k = 3; f(x, y) = x + y";
        f = function_3;
        break;
    case 4:
        f_name = "k = 4; f(x, y) = sqrt(x*x + y*y)";
        f = function_4;
        break;
    case 5:
        f_name = "k = 5; f(x, y) = x*x + y*y";
        f = function_5;
        break;
    case 6:
        f_name = "k = 6; f(x, y) = exp(x*x - y*y)";
        f = function_6;
        break;
    case 7:
        f_name = "k = 7; f(x, y) = 1 / (25(x*x + y*y) + 1)";
        f = function_7;
        break;
    }
    data.f = f;
}

// методы для рисования всего подряд
// в этих штуках меняются цвета
double search_k(double value, double abs_min, double abs_max) {
	if (value - abs_min < 0) return 0;
	if (abs_max - value < 0) return 1;
	if (fabs(abs_max - abs_min) < EPS) return 1;       // это для функции единица
	return fabs((value - abs_min) / (abs_max - abs_min));
}

void rgb_for_f(int &R, int &G, int &B, double k) {
	R = 127 + 127 * k;
	G = 127 - 127 * k;
	B = 127 - 127 * k;
}
void rgb_for_pf(int &R, int &G, int &B, double k) {
	R = 127 - 127 * k;
	G = 127 + 127 * k;
	B = 127 - 127 * k;
}
void rgb_for_residual(int &R, int &G, int &B, double k) {
	R = 127 - 127 * k;
	G = 127 - 127 * k;
	B = 127 + 127 * k;
}



QPointF Window::l2g(double x_loc, double y_loc, double y_min, double y_max) {
	double x_gl = (x_loc - data.a) / (data.b - data.a) * width();
	double y_gl = (y_max - y_loc)  / (y_max - y_min)   * height();
	return QPointF(x_gl, y_gl);
}

void Window::draw_one_element(QPainter *painter, double x1, double y1, double x2, double y2, double x3, double y3, QColor color) {
	QPainterPath path;
	//double min_y = data.c, max_y = data.d;

	path.moveTo(l2g(x1, y1, data.c, data.d));
	path.lineTo(l2g(x2, y2, data.c, data.d));
	path.lineTo(l2g(x3, y3, data.c, data.d));
	path.lineTo(l2g(x1, y1, data.c, data.d));

	painter->setPen(Qt::NoPen);
	painter->fillPath(path, QBrush(color));
}

void Window::draw_f(QPainter *painter) {
	double a = data.a, b = data.b, c = data.c, d = data.d;
	double mx = data.mx, my = data.my;
	double (*f)(double, double) = data.f;

	double hx = (b - a) / mx;
	double hy = (d - c) / my;
	//printf("F = %lf %lf. %lf %lf %lf %lf\n", hx, hy, a, b, c, d);
	double val_f, k_color;
	int R, G, B, i, j;
	
	double x1, y1, x2, y2, x3, y3;
	double abs_min = f(a + hx / 3.0, c + (2.0 / 3.0) * hy);
	double abs_max = abs_min;

	data.find_min_max(abs_min, abs_max);
	for (i = 0; i < mx; ++i) {
		for (j = 0; j < my; ++j) {
			x1 = a + i * hx;        y1 = c + j * hy;
			x2 = x1;                y2 = c + (j + 1) * hy;
			x3 = a + (i + 1) * hx;  y3 = y2;

			val_f = f(a + (i + 1.0 / 3.0) * hx, c + (j + 2.0 / 3.0) * hy);
			k_color = search_k(val_f, abs_min, abs_max);
			rgb_for_f(R, G, B, k_color);
			draw_one_element(painter, x1, y1, x2, y2, x3, y3, QColor(R, G, B));

			x2 = a + (i + 1) * hx;  y2 = c + j * hy;
			val_f = f(a + (i + 2.0 / 3.0) * hx, c + (j + 1.0 / 3.0) * hy);
			k_color = search_k(val_f, abs_min, abs_max);
			rgb_for_f(R, G, B, k_color);
			draw_one_element(painter, x1, y1, x2, y2, x3, y3, QColor(R, G, B));
		}
	}

	char max_str[100];
	double f_max = fmax(fabs(abs_min), fabs(abs_max));
	F_max = f_max;		printf("max{|F_min|, |F_max|} = %le\n", F_max);
	sprintf(max_str, "max{|F_min|, |F_max|} = %.2e", f_max);
	painter->setPen("black");
	painter->drawText(10, 40, max_str);
	draw_txt(painter);
}

void Window::draw_Pf(QPainter *painter) {
	double a = data.a, b = data.b, c = data.c, d = data.d;
	double mx = data.mx, my = data.my;

	double hx = (b - a) / mx;
	double hy = (d - c) / my;
	//printf("PF = %lf %lf. %lf %lf %lf %lf\n", hx, hy, a, b, c, d);
	double val_Pf, k_color;
	int R, G, B, i, j;
	
	double x1, y1, x2, y2, x3, y3;
	double abs_min = data.pf(a + hx / 3.0, c + (2.0 / 3.0) * hy);
	double abs_max = abs_min;

	data.pfind_min_max(abs_min, abs_max);
	for (i = 0; i < mx; ++i) {
		for (j = 0; j < my; ++j) {
			x1 = a + i * hx;        y1 = c + j * hy;
			x2 = x1;                y2 = c + (j + 1) * hy;
			x3 = a + (i + 1) * hx;  y3 = y2;

			val_Pf = data.pf(a + (i + 1.0 / 3.0) * hx, c + (j + 2.0 / 3.0) * hy);
			k_color = search_k(val_Pf, abs_min, abs_max);
			rgb_for_pf(R, G, B, k_color);
			draw_one_element(painter, x1, y1, x2, y2, x3, y3, QColor(R, G, B));

			x2 = a + (i + 1) * hx;  y2 = c + j * hy;
			val_Pf = data.pf(a + (i + 2.0 / 3.0) * hx, c + (j + 1.0 / 3.0) * hy);
			k_color = search_k(val_Pf, abs_min, abs_max);
			rgb_for_pf(R, G, B, k_color);
			draw_one_element(painter, x1, y1, x2, y2, x3, y3, QColor(R, G, B));
		}
	}

	char max_str[100];
	double pf_max = fmax(fabs(abs_min), fabs(abs_max));
	PF_max = pf_max;		printf("max{|F_min|, |F_max|} = %le\n", PF_max);
	sprintf(max_str, "max{|F_min|, |F_max|} = %.2e", pf_max);
	painter->setPen("black");
	painter->drawText(10, 40, max_str);
	draw_txt(painter);
}

void Window::draw_residual(QPainter *painter) {
	double a = data.a, b = data.b, c = data.c, d = data.d;
	double mx = data.mx; double my = data.my; double (*f)(double, double) = data.f;

	double hx = (b - a) / mx;
	double hy = (d - c) / my;
	double val_res, f_val, pf_val, k_color;
	int R, G, B, i, j;
	
	double x1, y1, x2, y2, x3, y3;
	double abs_min = fabs(f(a + hx / 3.0, c + (2.0 / 3.0) * hy) - data.pf(a + hx / 3.0, c + (2.0 / 3.0) * hy));
	double abs_max = abs_min;
	data.residual_min_max(abs_min, abs_max);

	for (i = 0; i < mx; ++i) {
		for (j = 0; j < my; ++j) {
			x1 = a + i * hx;        y1 = c + j * hy;
			x2 = x1;                y2 = c + (j + 1) * hy;
			x3 = a + (i + 1) * hx;  y3 = y2;

			f_val = f(a + (i + 1.0 / 3.0) * hx, c + (j + 2.0 / 3.0) * hy);
			pf_val = data.pf(a + (i + 1.0 / 3.0) * hx, c + (j + 2.0 / 3.0) * hy);

			val_res = fabs(f_val - pf_val);
			k_color = search_k(val_res, abs_min, abs_max);
			rgb_for_residual(R, G, B, k_color);
			draw_one_element(painter, x1, y1, x2, y2, x3, y3, QColor(R, G, B));
			
			x2 = a + (i + 1) * hx;	y2 = c + j * hy;
			f_val = f(a + (i + 2.0 / 3.0) * hx, c + (j + 1.0 / 3.0) * hy);
			pf_val = data.pf(a + (i + 2.0 / 3.0) * hx, c + (j + 1.0 / 3.0) * hy);

			val_res = fabs(f_val - pf_val);
			k_color = search_k(val_res, abs_min, abs_max);
			rgb_for_residual(R, G, B, k_color);
			draw_one_element(painter, x1, y1, x2, y2, x3, y3, QColor(R, G, B));
		}
	}

	char max_str[100];
	double res_max = fmax(fabs(abs_min), fabs(abs_max));
	RES_max = res_max;		printf("max{|F_min|, |F_max|} = %le\n", RES_max);
	sprintf(max_str, "max{|f - Pf|} = %.2e", res_max);
	painter->setPen("black");
	painter->drawText(10, 40, max_str);
	draw_txt(painter);
}



// ФУНКЦИИ ДЛЯ РАБОТЫ В ОКНЕ

void Window::change_func() {
	if (msr_ready()) {
		func_id = (func_id + 1) % 8;
		data.func_id = func_id;
		set_f();

		double abs_min = DEF_ABS_MIN, abs_max = DEF_ABS_MIN;
		f_max = data.find_min_max(abs_min, abs_max);
		//printf("f_max = %le\n", f_max);
		data.norma = /* (int) */f_max;

		data.repeat_new_data();

		for (int k = 0; k < data.p; k++) {
			args[k].copy_data(data);
			
			args[k].k = k;
			args[k].mutex = &mutex;
			args[k].cond = &cond;
			args[k].ready = false;
		}

		if (threads_created) pthread_cond_broadcast(&cond);
		else {
			threads_created = true;
			for (int i = 0; i < data.p; i++) pthread_create(&args[i].tid, 0, thread_func, args + i);
		}
		QTimer::singleShot(20, this, &Window::msr_wait);
        //one_deletion_of_results();
	}
	else QMessageBox::warning(0, "Внимание", "Подождите конца вычислений.");
}

void Window::change_mode() {
	mode = (mode + 1) % 3;
	update();
	// то что ниже работает плохо
	/* switch (mode) {
		case 0:		printf("max{|F_min|, |F_max|} = %le\n", F_max); break;
		case 1:		printf("max{|F_min|, |F_max|} = %le\n", PF_max); break;
		case 2:		printf("max{|F_min|, |F_max|} = %le\n", RES_max); break;
	} */
}

void Window::twice_n() {
	if (msr_ready()) {
		data.nx *= 2;
		data.ny *= 2;
		data.repeat_new_data();

		for (int i = 0; i < data.p; i++) {
			args[i].copy_data(data);
			args[i].k = i;
			args[i].mutex = &mutex;
			args[i].cond = &cond;
			args[i].ready = false;
		}
		pthread_cond_broadcast(&cond);
		QTimer::singleShot(20, this, &Window::msr_wait);
        //one_deletion_of_results();
	}
	else QMessageBox::warning(0, "Внимание", "Подождите конца вычислений.");
}

void Window::halve_n() {
	if (msr_ready()) {
		if (data.nx >= 4) data.nx /= 2;
		if (data.ny >= 4) data.ny /= 2;
		data.repeat_new_data();

		for (int i = 0; i < data.p; i++) {
			args[i].copy_data(data);
			args[i].k = i;
			args[i].mutex = &mutex;
			args[i].cond = &cond;
			args[i].ready = false;
		}

		pthread_cond_broadcast(&cond);
		QTimer::singleShot(20, this, &Window::msr_wait);
        //one_deletion_of_results();
	} 
	else QMessageBox::warning(0, "Внимание", "Подождите конца вычислений.");
}

void Window::increase_p() {
	if (msr_ready()) {
		parameter += 1;
		data.repeat_new_data();
		data.parameter += 1;

		for (int i = 0; i < data.p; i++) {
			args[i].copy_data(data);
			args[i].k = i;
			args[i].mutex = &mutex;
			args[i].cond = &cond;
			args[i].ready = false;
		}

		pthread_cond_broadcast(&cond);
		QTimer::singleShot(20, this, &Window::msr_wait);
        //one_deletion_of_results();
	}
	else QMessageBox::warning(0, "Внимание", "Подождите конца вычислений.");
}

void Window::decrease_p() {
	if (msr_ready()) {
		parameter -= 1;
		data.repeat_new_data();
		data.parameter -= 1;

		for (int i = 0; i < data.p; i++) {
			args[i].copy_data(data);
			args[i].k = i;
			args[i].mutex = &mutex;
			args[i].cond = &cond;
			args[i].ready = false;
		}

		pthread_cond_broadcast(&cond);
		QTimer::singleShot(20, this, &Window::msr_wait);
        //one_deletion_of_results();
	}
	else QMessageBox::warning(0, "Внимание", "Подождите конца вычислений.");
}

void Window::zoom_in() {
	double len_x = data.b - data.a, len_y = data.d - data.c;

	parameter = 0;
	data.parameter = 0;

	data.a += len_x / 4;
	data.b -= len_x / 4;
	data.c += len_y / 4;
	data.d -= len_y / 4;

	//halve_m();
	update();
}

void Window::zoom_out() {
	double len_x = data.b - data.a, len_y = data.d - data.c;

	parameter = 0;
	data.parameter = 0;

	data.a -= len_x / 2;
	data.b += len_x / 2;
	data.c -= len_y / 2;
	data.d += len_y / 2;

	//twice_m();
	update();
}

void Window::twice_m() {
	data._mx = data._mx * 2;	data._my = data._my * 2;
	/* if (data.mx <= 320) data.mx *= 2;
	if (data.my <= 180) data.my *= 2; */
	data.mx *= 2;	data.my *= 2;
	if (data.mx > DEF__MX)	data.mx = DEF__MX;
	if (data.my > DEF__MY)	data.my = DEF__MY;

	double abs_min = DEF_ABS_MIN, abs_max = DEF_ABS_MIN;
	f_max = data.find_min_max(abs_min, abs_max);

	update();
}

void Window::halve_m() {
	// если _mx и _my деленный на два больше DEF__MX и DEF__MY то его считаем таким
	if (data._mx >= 4) data._mx = data._mx / 2;
	if (data._my >= 4) data._my = data._my / 2;
	if (data._mx > DEF__MX) data.mx = DEF__MX;
	else if (data.mx >= 4)	data.mx = data._mx;
	if (data._my > DEF__MY) data.my = DEF__MY;
	else if (data.my >= 4)	data.my = data._my;

	double abs_min = DEF_ABS_MIN, abs_max = DEF_ABS_MIN;
	f_max = data.find_min_max(abs_min, abs_max);

	update();
}



// ЕЩЕ ДВА МЕТОДА

void Window::close() {
	if (msr_ready()) {
		//one_deletion_of_results();
		/* if (args)	delete[] args;
		if (tid)	delete[] tid; */
		widget->close();
	}
	else QMessageBox::warning(0, "Внимание", "Подождите конца вычислений.");
}

void Window::paintEvent(QPaintEvent *) {
	QPainter painter(this);
	if (mode == 0)       draw_f(&painter);
	else if (mode == 1)  draw_Pf(&painter);
	else                 draw_residual(&painter);
}

