#include "glut.h"
#include <GL/gl.h>
#include <vector>
#include <array>
#include <cmath>

constexpr float blue[] = { 0.0, 0.0, 1.0, 1.0 };
constexpr float white[] = { 1.0, 1.0, 1.0, 1.0 };
constexpr float white_and_transparent[] = { 1.0, 1.0, 1.0, 0.5 };
constexpr float black[] = { 0.0, 0.0, 0.0, 1.0 };
constexpr float transparent[] = { 0.0, 0.0, 0.0, 0.0 };

class MassPoint
{
public:
    double height;
    double velocity;
public:
    MassPoint(double y = 0, double v = 0)
    {
		height = y;
		velocity = v;
    }

    void addVelocity(double v)
    {
		velocity += v;
	}

    void update(double dv, double dt)
    {
        velocity += dv - velocity * 0.001;
        velocity -= height * 0.0001;
//        velocity += dv;
        height += velocity * dt;
    }
};

struct Params {
    float p[4];
    Params() = default;
    Params(float const* p) {
        for (int i = 0; i < 4; i++)
        {
            this->p[i] = p[i];
        }
    }

    Params& maltiply(float a) {
        for (int i = 0; i < 3; i++)
        {
			p[i] *= a;
		}
		return *this;
	}

    Params& alpha(float a) {
        p[3] = a;
        return *this;
    }
    operator float* () { return p; }
};

struct Material {
	Params ambient;
	Params diffuse;
	Params specular;
	float shininess;
	Material() = default;
	Material(Params const& a, Params const& d, Params const& s, float sh) : ambient(a), diffuse(d), specular(s), shininess(sh) {}
    void apply() {
		glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
		glMaterialf(GL_FRONT, GL_SHININESS, shininess);
	}
};

class WaterSurface
{
    std::vector<std::vector<MassPoint*>> points;
    double tension;
public:
    Material material = {};

    WaterSurface(uint64_t n, uint64_t m, double t)
    {
		tension = t;
		points.resize(n);
        for (uint64_t i = 0; i < n; i++)
        {
			points[i].resize(m);
            for (uint64_t j = 0; j < m; j++)
            {
				points[i][j] = new MassPoint();
			}
		}
	}

    void addForce(uint64_t i, uint64_t j, double f)
    {
		points[i][j]->addVelocity(f);
    }

    void setHeight(uint64_t i, uint64_t j, double h)
    {
        points[i][j]->height = h;
    }

    void update(double dt)
    {
        std::vector<std::vector<double>> force(points.size(), std::vector<double>(points[0].size(), 0));

        for (size_t i = 0; i < points.size() - 1; i++) {
            for (size_t j = 0; j < points[0].size(); j++) {
                force[i][j] += points[i + 1][j]->height;
                force[i+1][j] += points[i][j]->height;

                force[i][j] -= points[i][j]->height;
                force[i+1][j] -= points[i+1][j]->height;
            }
        }

        for (size_t i = 0; i < points.size(); i++) {
            for (size_t j = 0; j < points[0].size() - 1; j++) {
                force[i][j] += points[i][j + 1]->height;
                force[i][j + 1] += points[i][j]->height;

                force[i][j] -= points[i][j]->height;
                force[i][j + 1] -= points[i][j + 1]->height;
            }
        }

        for (size_t i = 0; i < points.size(); i++) {
            for (size_t j = 0; j < points[0].size(); j++) {
				points[i][j]->update(force[i][j] * tension, dt);
			}
        }
    }

    static constexpr double NORMAL_Y = 1;

    void draw()
    {

        std::vector<std::vector<std::array<double, 2>>> normals(points.size(), std::vector<std::array<double, 2>>(points[0].size(), {0, 0}));
        for (size_t i = 1; i < points.size() - 1; i++) {
            for (size_t j = 1; j < points[0].size() - 1; j++) {
                normals[i][j][0] = points[i - 1][j]->height - points[i + 1][j]->height;
                normals[i][j][1] = points[i][j - 1]->height - points[i][j + 1]->height;
            }
        }

        size_t last_x = points.size() - 1;
        size_t last_y = points[0].size() - 1;

        for (size_t i = 1; i < points.size() - 1; i++) {
            normals[i][0][0] = points[i - 1][0]->height - points[i + 1][0]->height;
            normals[i][0][1] = (points[i][0]->height - points[i][1]->height) * 2;

            normals[i][last_y][0] = points[i - 1][last_y]->height - points[i + 1][last_y]->height;
            normals[i][last_y][1] = (points[i][last_y - 1]->height - points[i][last_y]->height) * 2;
        }

        for (size_t i = 1; i < points[0].size() - 1; i++) {
			normals[0][i][0] = (points[0][i]->height - points[1][i]->height) * 2;
			normals[0][i][1] = points[0][i - 1]->height - points[0][i + 1]->height;

			normals[last_x][i][0] = (points[last_x - 1][i]->height - points[last_x][i]->height) * 2;
			normals[last_x][i][1] = points[last_x][i - 1]->height - points[last_x][i + 1]->height;
        }

        normals[0][0][0] = (points[0][0]->height - points[1][0]->height) * 2;
        normals[0][0][1] = (points[0][0]->height - points[0][1]->height) * 2;

        normals[0][last_y][0] = (points[0][last_y]->height - points[1][last_y]->height) * 2;
        normals[0][last_y][1] = (points[0][last_y - 1]->height - points[0][last_y]->height) * 2;

        normals[last_x][0][0] = (points[last_x - 1][0]->height - points[last_x][0]->height) * 2;
        normals[last_x][0][1] = (points[last_x][0]->height - points[last_x][1]->height) * 2;

        normals[last_x][last_y][0] = (points[last_x - 1][last_y]->height - points[last_x][last_y]->height) * 2;
        normals[last_x][last_y][1] = (points[last_x][last_y - 1]->height - points[last_x][last_y]->height) * 2;

        double center_x = (points.size() - 1) / 2.0;
        double center_z = (points[0].size() - 1) / 2.0;

        material.apply();

        for (size_t i = 0; i < points.size() - 1; i++) {
            glBegin(GL_TRIANGLE_STRIP);
            for (size_t j = 0; j < points[0].size(); j++) {
                glNormal3d(normals[i][j][0], NORMAL_Y, normals[i][j][1]);
				glVertex3d((i - center_x) * 0.1, points[i][j]->height, (j - center_z) * 0.1);

				glNormal3d(normals[i + 1][j][0], NORMAL_Y, normals[i + 1][j][1]);
				glVertex3d((i + 1 - center_x) * 0.1, points[i + 1][j]->height, (j - center_z) * 0.1);
            }
            glEnd();
		}
    }
};

void display();
void update(int timer_id);
void drawWaterSurface(int timer_id);
void rain(int timer_id);

WaterSurface* waterSurface;

int main(int argc, char** argv)
{
    glutInit(&argc, argv);

    waterSurface = new WaterSurface(100, 100, 0.1);
    waterSurface->material = Material(black, Params(white).alpha(0.1), Params(white).maltiply(0.000000000000000000001f).alpha(0.5), 128);

    glutInitWindowPosition(0, 0);
    glutInitWindowSize(800, 800);
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

    glTranslated(0, 0, -5);

    // ”wŒi
    glMaterialfv(GL_FRONT, GL_DIFFUSE, Params(white).maltiply(0.5));
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

    glTranslated(0, 0, -5);

    if (waterSurface != nullptr)
    {
        glPushMatrix();
        glScaled(1, 1, 1);
        waterSurface->draw();
        glPopMatrix();
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