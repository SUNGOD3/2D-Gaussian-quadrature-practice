//#define _CRT_SECURE_NO_WARNINGS //dev c++
#include <iostream>
#include <iomanip>
#include <cmath>
//#include <GL/gl.h> //dev c++
#include <GL/glut.h>


//max number of interval
#define N 4
//max number of Gaussian points
#define K 4
using namespace std;
//Sample points in [-1, 1], 1-, 2-, 3-, and 4-points. Ignore extra points (the 0s).
long double  p[K][K]={{0.0, 0.0, 0.0, 0.0}, {0.5773502691896257, -0.5773502691896257, 0.0, 0.0}, {0.0, 0.7745966692414834, -0.7745966692414834, 0.0}, {0.3399810435848563, -0.3399810435848563, 0.8611363115940526, -0.8611363115940526}};
//function values in N interval
long double  Ji[K]={1,0.25,0.111111111111111111,0.0625},ans[K][K];
long double  wgt[K][K]={{2.0, 0.0, 0.0, 0.0}, {1.0, 1.0, 0.0, 0.0}, {0.888888888888888888888, 0.555555555555555556, 0.555555555555555556, 0.0},{0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538}};
//integral range = [0, 2].
long double   a=0.0, b=2.0,pi=3.14159265358979323846;

long double my_func(long double x, long double y){
	return sin(2.0*pi*x)*cos(3.0*pi*y)+1.0;
}

long double conform_map(long double t, long double a, long double b){
   return (b-a)*t/2.0 + (a+b)/2.0;
}

long double  gauss_2D_quadrature(int k,long double ax,long double bx,long double ay,long double by){
    long double   sum=0.0 , x , y;
    for(int i=0;i<k;++i){
    	for(int j=0;j<k;++j){
        	x = conform_map(p[k-1][i], ax, bx); //Conformation mapping
        	y = conform_map(p[k-1][j], ay, by);
			sum += wgt[k-1][i]*wgt[k-1][j]*my_func(x,y); //The integration
    	}
	}
    //sum = sum*(b-a)/2.0; //Scale the integration value.
    //sum = sum*(b-a)/2.0;
    return sum;
}
void draw_unit(int k) {
	//畫x,y軸
	glBegin(GL_LINES);
	glColor3f(0.5, 0.5, 0.5);
	glVertex2f(-1.0, 0.0);
	glVertex2f(1.0, 0.0);
	glVertex2f(0.0, 1.0);
	glVertex2f(0.0, -1.0);
	glEnd();
	//畫格子
	glBegin(GL_LINES);
	glColor3f(0.8, 0.8, 0.8);
	for (float i = -20; i < 20; ++i) {
		glVertex2f(-1.0, i / 20);
		glVertex2f(1.0, i / 20);
		glVertex2f(i / 20, 1.0);
		glVertex2f(i / 20, -1.0);
	}
	glEnd();
	//點塗色
	glPointSize(5);    //點的大小
	glBegin(GL_POINTS); //選擇畫點
	for (int i = 0; i < N; ++i) {
			glColor3f(1, 0, 0);
			glVertex2f((long double)i / 10.0 , ans[k-1][i] * 10000000000000000 / 15.0);
	}
	glEnd();
	//畫單位點
	glPointSize(3);    //點的大小
	glBegin(GL_POINTS); //選擇畫點
	glColor3f(0, 0, 0);
	for (float i = -20; i < 20; i++) {
		glVertex2f(0.0, i / 20);
		glVertex2f(i / 20, 0.0);
	}
	glEnd();

}
void display_func() {
	glClearColor(1, 1, 1, 1);//the color of eraser
	glClear(GL_COLOR_BUFFER_BIT);//擦乾淨
	glViewport(0, 350, 350, 350);//左上
	
	draw_unit(1);
	glViewport(350, 350, 350, 350);//右上
	draw_unit(2);
	glViewport(0,0,350,350);//左下
	draw_unit(3);
	glViewport(350, 0, 350, 350);//右下
	draw_unit(4);
	glFlush();
}

void reshape(int w, int h) {
	glViewport(0, 0, w, h);    //調整視角
}

int main(int argc, char** argv) {
  int        k, numInterval, i, j;
  long double  gaussSum, localSum,exactSum, h;
  glutInit(&argc, argv);
  glutInitWindowPosition(500, 0);    //初始位置
  glutInitWindowSize(700, 700);    //大小
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);    //設定模式
  glutCreateWindow("Test1");
  glutReshapeFunc(reshape);
  glutDisplayFunc(display_func);
  	cout << "Q2:\n"; 
  //Using different number of Gaussian points
  	for(k=1;k<=4;k++){
		cout << "No. of Gaussian points = " << k << "\nNo. of ceil = ";
		for(int i=1;i<=4;++i){
			cout << i*i << setw(23) ;
		}
		cout << "\n"; 
	  	//Using different numbers of intervals
	    for(numInterval=1; numInterval<=N; numInterval++){
		    //printf("  No. of interval = %d\n", numInterval);
	        h = (b-a)/numInterval;
		    gaussSum = 0.0;
		    for(long double  i=0;i<numInterval;++i){
		      	for(long double j=0;j<numInterval;++j){
	           		localSum = gauss_2D_quadrature(k, a+i*h, a+(i+1)*h , a+j*h, a+(j+1)*h);
			   		gaussSum += localSum;
			  	}
			}
			gaussSum*=Ji[numInterval-1];
			cout << setw(12) << fixed <<  setprecision(20) << gaussSum << " ";
			ans[k-1][numInterval-1]=fabs(4.0-gaussSum);
		}
		cout << "\nRelative Error\n";
		for(int i=0;i<N;++i){
			cout << setw(12) << fixed  <<  setprecision(20) << ans[k-1][i] << " ";
		}
		cout << "\n";
  	}
  	cout << "\nQ4:\n隨著cell遞增，可以看出除了cell = 4的情況外，基本上準確度是隨之遞增的。";
  	cout << "\nQ5:\n隨著取樣點個數增加，準確度的變化似乎不太明顯。";
  	cout << "\nQ6:\ncell";
	glutMainLoop();
} 

