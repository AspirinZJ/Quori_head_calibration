#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h> // using glut
// #include "rt_nonfinite.h"
#include "sphereToimage.h"
#include <stddef.h>
#include "rotation_theta.h"

const int SCREENWIDTH = 1920.0;
const int SCREENHEIGHT = 1080.0;
// const double RADIUS = 3;  // the radius of the sphere screen
const double RADIUS = 4.0;          // the radius of the sphere screen
double dx = -2.51;                  // the parameter for calibrating the x-axis installation error
double dy = -1.14;                  // the parameter for calibrating the y-axis installation error
int Glut_Handle1_Window_FullScreen; // Glut create handle

void displayFunc();
void idleFunc();
void kbdFunc(unsigned char key, int x, int y);

int main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // double buffering
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(SCREENWIDTH, SCREENHEIGHT);
    Glut_Handle1_Window_FullScreen = glutCreateWindow("Full Screen Window");

    glClearColor(0.0, 0.0, 0.0, 0.0); // black background
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    glutFullScreen(); // full screen mode
    glutDisplayFunc(displayFunc);
    glutIdleFunc(idleFunc);
    glutKeyboardFunc(kbdFunc);
    glutMainLoop();

    return 0;
}

void displayFunc() // display function
{
    float Pi = 3.141592653;
    // double theta0 = (double)1 / 2 * Pi;
    double theta0 = (double)0.95 / 2 * Pi; // initial polar angle, on the equator
    double psi0;                           // initial azimuthal angle
    double dtheta = theta0 - Pi / 2;
    double px, py;
    double semiMajor = 0.3, semiMinor = 0.3; // semi-major and semi-minor
    double phi = 0;
    int polyNum = 100; // use a 100-edge polygon to approximate a circle
    double CoordX[100], CoordY[100];
    double dPi = 2 * Pi / polyNum; // Pi step size
    double x, y, xRot, yRot;       // x and y (and rotated xRot and yRot)position of a circle
    double theta, psi, Theta, Psi;
    double t; // parametric equation of a circle

    glColor3f(1.0, 0.2, 0.2); // color of the circle
    glClear(GL_COLOR_BUFFER_BIT);

    // draw the circles
    for (float theta0 = 0.22 * Pi; theta0 < 0.9 * Pi; theta0 += 0.12 * Pi)
    {
        for (float psi0 = 0; psi0 < 2 * Pi; psi0 += 0.5 * Pi)
        {
            for (int i = 0; i < polyNum; i++)
            {
                t = i * dPi;
                // x and y save the position of an ellipse
                x = semiMajor * cos(t);
                y = semiMinor * sin(t);
                // x and y time a rotation matrix
                xRot = x * cos(phi) - y * sin(phi);
                yRot = x * sin(phi) + y * cos(phi);
                // approximate an circle from a flat plane to a sphere
                psi = xRot / (RADIUS * sin(theta0)) + psi0;
                theta = theta0 - yRot / RADIUS;
                rotation_theta(theta, psi, psi0, dtheta, &Theta, &Psi);                 // rotate the circle from the equator to other places
                sphereToimage(Theta, Psi, SCREENHEIGHT, SCREENWIDTH, &px, &py, dx, dy); // warp the 3d spherical position onto a 2D image plane
                px = px - SCREENWIDTH / 2;                                              // pixel coordinates
                py = py - SCREENHEIGHT / 2;                                             // pixel coordinates
                // py = -py;
                CoordX[i] = 2 * px / SCREENHEIGHT;
                CoordY[i] = 2 * py / SCREENHEIGHT;
            }
            glBegin(GL_POLYGON);
            for (int i = 0; i < polyNum; i++)
                glVertex2f(CoordX[i] * SCREENHEIGHT / SCREENWIDTH, CoordY[i]);
            glEnd();
        }
    }

    glFlush();
    glutSwapBuffers(); // for double buffering
}

void idleFunc() // idle function
{
    glutPostRedisplay(); // rerender(call displayFunc())
}

void kbdFunc(unsigned char key, int x, int y)
{
    if ('w' == key || 'W' == key)
    {
        dy += 0.01;
        printf("dx = %lf, dy = %lf\n", dx, dy);
    }
    if ('s' == key || 'S' == key)
    {
        dy -= 0.01;
        printf("dx = %lf, dy = %lf\n", dx, dy);
    }
    if ('a' == key || 'A' == key)
    {
        dx -= 0.01;
        printf("dx = %lf, dy = %lf\n", dx, dy);
    }
    if ('d' == key || 'D' == key)
    {
        dx += 0.01;
        printf("dx = %lf, dy = %lf\n", dx, dy);
    }

    if ('q' == key || 'Q' == key)
    {
        printf("Quitting the program.\n");
        glutDestroyWindow(Glut_Handle1_Window_FullScreen);
        exit(0);
    }
}