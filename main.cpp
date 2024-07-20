#include "glut.h"
#include <GL/gl.h>
#include <vector>
#include <array>
#include <cmath>

float blue[] = { 0.0, 0.0, 1.0, 1.0 };
float white[] = { 1.0, 1.0, 1.0, 1.0 };
float white_and_transparent[] = { 1.0, 1.0, 1.0, 0.5 };
float black[] = { 0.0, 0.0, 0.0, 1.0 };
float transparent[] = { 0.0, 0.0, 0.0, 0.0 };

class MassPoint
{
public:
    double mass;
    double position[3];
    double velocity;
public:
    MassPoint(double m, double x, double y, double z, double v = 0)
    {
		mass = m;
		position[0] = x;
		position[1] = y;
		position[2] = z;
		velocity = v;
    }

    void addVelocity(double v)
    {
		velocity += v;
	}

    void update(double dv, double dt)
    {
        velocity += dv - velocity * 0.001;
        velocity -= position[1] * 0.0001;
//        velocity += dv;
        position[1] += velocity * dt;
    }
};

class WaterSurface
{
    std::vector<std::vector<MassPoint*>> points;
    double tension;
public:
    WaterSurface(uint64_t n, uint64_t m, double t)
    {
		tension = t;
		points.resize(n);
        for (uint64_t i = 0; i < n; i++)
        {
			points[i].resize(m);
            for (uint64_t j = 0; j < m; j++)
            {
				points[i][j] = new MassPoint(1.0, i * 0.1, 0, j * 0.1);
			}
		}
	}

    void addForce(uint64_t i, uint64_t j, double f)
    {
		points[i][j]->addVelocity(f);
    }

    void setHeight(uint64_t i, uint64_t j, double h)
    {
        points[i][j]->position[1] = h;
    }

    void update(double dt)
    {
        std::vector<std::vector<double>> force(points.size(), std::vector<double>(points[0].size(), 0));

        for (size_t i = 0; i < points.size() - 1; i++) {
            for (size_t j = 0; j < points[0].size(); j++) {
                force[i][j] += points[i + 1][j]->position[1];
                force[i+1][j] += points[i][j]->position[1];

                force[i][j] -= points[i][j]->position[1];
                force[i+1][j] -= points[i+1][j]->position[1];
            }
        }

        for (size_t i = 0; i < points.size(); i++) {
            for (size_t j = 0; j < points[0].size() - 1; j++) {
                force[i][j] += points[i][j + 1]->position[1];
                force[i][j + 1] += points[i][j]->position[1];

                force[i][j] -= points[i][j]->position[1];
                force[i][j + 1] -= points[i][j + 1]->position[1];
            }
        }

        for (size_t i = 0; i < points.size(); i++) {
            for (size_t j = 0; j < points[0].size(); j++) {
				points[i][j]->update(force[i][j] * tension, dt);
			}
        }
    }

    void draw()
    {

        std::vector<std::vector<std::array<double, 3>>> normals(points.size(), std::vector<std::array<double, 3>>(points[0].size(), {0, 0, 0}));
        for (size_t i = 1; i < points.size() - 1; i++) {
            for (size_t j = 1; j < points[0].size() - 1; j++) {
                normals[i][j][0] = points[i - 1][j]->position[1] - points[i + 1][j]->position[1];
                normals[i][j][1] = 2;
                normals[i][j][2] = points[i][j - 1]->position[1] - points[i][j + 1]->position[1];
            }
        }

        for (size_t i = 1; i < points.size() - 1; i++) {
            normals[i][0][0] = points[i - 1][0]->position[1] - points[i + 1][0]->position[1];
            normals[i][0][1] = 2;
            normals[i][0][2] = (points[i][0]->position[1] - points[i][1]->position[1]) * 2;

            size_t last = points[0].size() - 1;
            normals[i][last][0] = points[i - 1][last]->position[1] - points[i + 1][last]->position[1];
            normals[i][last][1] = 2;
            normals[i][last][2] = (points[i][last - 1]->position[1] - points[i][last]->position[1]) * 2;
        }

        for (size_t i = 1; i < points[0].size() - 1; i++) {
			normals[0][i][0] = (points[0][i]->position[1] - points[1][i]->position[1]) * 2;
			normals[0][i][1] = 2;
			normals[0][i][2] = points[0][i - 1]->position[1] - points[0][i + 1]->position[1];

			size_t last = points.size() - 1;
			normals[last][i][0] = (points[last][i]->position[1] - points[last - 1][i]->position[1]) * 2;
			normals[last][i][1] = 2;
			normals[last][i][2] = points[last][i - 1]->position[1] - points[last][i + 1]->position[1];
        }

        normals[0][0][0] = (points[0][0]->position[1] - points[1][0]->position[1]) * 2;
        normals[0][0][1] = 2;
        normals[0][0][2] = (points[0][0]->position[1] - points[0][1]->position[1]) * 2;

        normals[0][points[0].size() - 1][0] = (points[0][points[0].size() - 1]->position[1] - points[1][points[0].size() - 1]->position[1]) * 2;
        normals[0][points[0].size() - 1][1] = 2;
        normals[0][points[0].size() - 1][2] = (points[0][points[0].size() - 2]->position[1] - points[0][points[0].size() - 1]->position[1]) * 2;

        normals[points.size() - 1][0][0] = (points[points.size() - 1][0]->position[1] - points[points.size() - 2][0]->position[1]) * 2;
        normals[points.size() - 1][0][1] = 2;
        normals[points.size() - 1][0][2] = (points[points.size() - 1][1]->position[1] - points[points.size() - 1][0]->position[1]) * 2;

        normals[points.size() - 1][points[0].size() - 1][0] = (points[points.size() - 1][points[0].size() - 1]->position[1] - points[points.size() - 2][points[0].size() - 1]->position[1]) * 2;
        normals[points.size() - 1][points[0].size() - 1][1] = 2;
        normals[points.size() - 1][points[0].size() - 1][2] = (points[points.size() - 1][points[0].size() - 2]->position[1] - points[points.size() - 1][points[0].size() - 1]->position[1]) * 2;

        
        for (size_t i = 0; i < points.size() - 1; i++) {
            glBegin(GL_TRIANGLE_STRIP);
            for (size_t j = 0; j < points[0].size(); j++) {
                glNormal3d(normals[i][j][0], normals[i][j][1], normals[i][j][2]);
				glVertex3d(points[i][j]->position[0], points[i][j]->position[1], points[i][j]->position[2]);

				glNormal3d(normals[i + 1][j][0], normals[i + 1][j][1], normals[i + 1][j][2]);
				glVertex3d(points[i + 1][j]->position[0], points[i + 1][j]->position[1], points[i + 1][j]->position[2]);
            }
            glEnd();
		}
    }
};

WaterSurface* waterSurface;
struct Params {
    float p[4];
    Params() = default;
    Params(float* p) {
        for (int i = 0; i < 4; i++)
        {
            this->p[i] = p[i];
        }
    }

    Params& alpha(float a) {
		p[3] = a;
		return *this;
	}
    operator float*() { return p; }
};

Params multiply(Params a, double value)
{
	Params result;
    result.p[3] = a.p[3];
    for (int i = 0; i < 3; i++)
    {
		result.p[i] = a.p[i] * value;
	}
	return result;

}
void display();
void update(int timer_id);
void drawWaterSurface(int timer_id);
void rain(int timer_id);

int main(int argc, char** argv)
{
    glutInit(&argc, argv);

    waterSurface = new WaterSurface(100, 100, 0.1);
    
    waterSurface->setHeight(50, 50, 2.0);


    glutInitWindowPosition(0, 0);
    glutInitWindowSize(400, 400);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_ALPHA);
    glutCreateWindow("“d20-0070 ¬–x³Ž÷");
    glClearColor(0.0, 0.0, 0.0, 1.0);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, 1.0, 1.0, 50.0);
//    gluLookAt(0, 1, 0, 0, 0.5, -8, 0, 1, 0);
    gluLookAt(0, 10, 0, 0, 0, -8, 0, 1, 0);

    glShadeModel(GL_SMOOTH);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);

    float light_position[] = { 1.0, 1.0, -1.0, 0.0};
    float white_and_transparent[] = { 1.0, 1.0, 1.0, 0.5 };
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
    glLightfv(GL_LIGHT0, GL_SPECULAR, white);
    glLightfv(GL_LIGHT0, GL_AMBIENT, white);

    glutDisplayFunc(display);
    glutTimerFunc(100, update, 0);
    glutTimerFunc(100, drawWaterSurface, 1);
    glutTimerFunc(100, rain, 2);
    glutMainLoop();
	return 0;
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

    // ”wŒi
    glMaterialfv(GL_FRONT, GL_DIFFUSE, multiply(white, 0.5));
    glMaterialfv(GL_FRONT, GL_SPECULAR, black);
    glMaterialfv(GL_FRONT, GL_AMBIENT, black);
    glMaterialf(GL_FRONT, GL_SHININESS, 0);
    glBegin(GL_QUADS);
    glNormal3d(0, 1, 0);
    glVertex3d(-10, -1, -10);
    glVertex3d(-10, -1, 10);
    glVertex3d(10, -1, 10);
    glVertex3d(10, -1, -10);
    glEnd();

    glTranslated(-5, 0, -10);

    if (waterSurface != nullptr)
    {
        glMaterialfv(GL_FRONT, GL_DIFFUSE, multiply(white_and_transparent, 1).alpha(0.1));
        glMaterialfv(GL_FRONT, GL_SPECULAR, multiply(white_and_transparent, 0.000000000000000000001f));
        //glMaterialfv(GL_FRONT, GL_SPECULAR, multiply(white_and_transparent, 0.00001f));
        glMaterialf(GL_FRONT, GL_SHININESS, 128);
        waterSurface->draw();
    }

    glFlush();

}

void drawWaterSurface(int timer_id)
{
    if (timer_id != 1)
    {
        return;
    }

    glutPostRedisplay();
    glutTimerFunc(50, drawWaterSurface, 1);
}

void update(int timer_id)
{
    if (timer_id != 0)
    {
		return;
	}
    if (waterSurface != nullptr)
    {
		waterSurface->update(0.01);
	}

	glutTimerFunc(1, update, 0);
}

void rain(int timer_id)
{
    if (timer_id != 2)
    {
		return;
	}

	waterSurface->addForce(rand() % 100, rand() % 100, -3);

	glutTimerFunc(50 + rand() % 100, rain, 2);
}