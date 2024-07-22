#include "glut.h"
#include <GL/gl.h>
#include <vector>
#include <array>
#include <cmath>
#include <iostream>

constexpr float FRAMES_PER_SECOND = 20;

constexpr float blue[] = { 0.0, 0.0, 1.0, 1.0 };
constexpr float white[] = { 1.0, 1.0, 1.0, 1.0 };
constexpr float white_and_transparent[] = { 1.0, 1.0, 1.0, 0.5 };
constexpr float black[] = { 0.0, 0.0, 0.0, 1.0 };
constexpr float transparent[] = { 0.0, 0.0, 0.0, 0.0 };

std::array<double, 3> light_position = { 0, 1, -30 };

class MassPoint
{
public:
    double height;
    double velocity;
    static constexpr double RESISTANCE = 0.005;
    static constexpr double GRAVITY = 0.001;
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
        velocity += dv - velocity * RESISTANCE;
        velocity -= height * GRAVITY;
//        velocity += dv;
        height += velocity * dt;
    }
};

struct Params {
    float p[4];
    Params() = default;
    Params(float r, float g, float b, float a) {
		p[0] = r;
		p[1] = g;
		p[2] = b;
		p[3] = a;
	}
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

inline float schlick(double c_spec, std::array<double, 3> const& L, std::array<double, 3> const& N) {
	double c = L[0] * N[0] + L[1] * N[1] + L[2] * N[2];

	return c_spec + (1 - c_spec) * pow(1 - c, 5);
}

class WaterSurface
{
    std::vector<std::vector<MassPoint*>> points;
    double tension;
    double specular;
public:
    Material material = {};

    WaterSurface(uint64_t n, uint64_t m, double t, double s = 0.02) : tension(t), specular(s)
    {
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

    static constexpr double NORMAL_Y = 2;

    void draw()
    {

        std::vector<std::vector<std::array<double, 3>>> normals(points.size(), std::vector<std::array<double, 3>>(points[0].size(), {0,NORMAL_Y, 0}));
        
        for (size_t i = 1; i < points.size() - 1; i++) {
            for (size_t j = 1; j < points[0].size() - 1; j++) {
                normals[i][j][0] = points[i - 1][j]->height - points[i + 1][j]->height;
                normals[i][j][2] = points[i][j - 1]->height - points[i][j + 1]->height;
            }
        }

        size_t last_x = points.size() - 1;
        size_t last_y = points[0].size() - 1;

        for (size_t i = 1; i < points.size() - 1; i++) {
            normals[i][0][0] = points[i - 1][0]->height - points[i + 1][0]->height;
            normals[i][0][2] = (points[i][0]->height - points[i][1]->height) * 2;

            normals[i][last_y][0] = points[i - 1][last_y]->height - points[i + 1][last_y]->height;
            normals[i][last_y][2] = (points[i][last_y - 1]->height - points[i][last_y]->height) * 2;
        }

        for (size_t i = 1; i < points[0].size() - 1; i++) {
			normals[0][i][0] = (points[0][i]->height - points[1][i]->height) * 2;
			normals[0][i][2] = points[0][i - 1]->height - points[0][i + 1]->height;

			normals[last_x][i][0] = (points[last_x - 1][i]->height - points[last_x][i]->height) * 2;
			normals[last_x][i][2] = points[last_x][i - 1]->height - points[last_x][i + 1]->height;
        }

        normals[0][0][0] = (points[0][0]->height - points[1][0]->height) * 2;
        normals[0][0][2] = (points[0][0]->height - points[0][1]->height) * 2;

        normals[0][last_y][0] = (points[0][last_y]->height - points[1][last_y]->height) * 2;
        normals[0][last_y][2] = (points[0][last_y - 1]->height - points[0][last_y]->height) * 2;

        normals[last_x][0][0] = (points[last_x - 1][0]->height - points[last_x][0]->height) * 2;
        normals[last_x][0][2] = (points[last_x][0]->height - points[last_x][1]->height) * 2;

        normals[last_x][last_y][0] = (points[last_x - 1][last_y]->height - points[last_x][last_y]->height) * 2;
        normals[last_x][last_y][2] = (points[last_x][last_y - 1]->height - points[last_x][last_y]->height) * 2;

        // ê≥ãKâª

        for (size_t i = 0; i < points.size(); i++) {
            for (size_t j = 0; j < points[0].size(); j++) {
                double len = sqrt(normals[i][j][0] * normals[i][j][0] + normals[i][j][1] * normals[i][j][1] + normals[i][j][2] * normals[i][j][2]);
				normals[i][j][0] /= len;
                normals[i][j][1] /= len;
				normals[i][j][2] /= len;
			}
		}
        

        double center_x = (points.size() - 1) / 2.0;
        double center_z = (points[0].size() - 1) / 2.0;

        material.apply();

        float buf[4] = { 0, 0, 0, 1 };

        //glEnable(GL_COLOR_MATERIAL);
        //glColorMaterial(GL_FRONT, GL_SPECULAR);

        for (size_t i = 0; i < points.size() - 1; i++) {
            glBegin(GL_TRIANGLE_STRIP);
            for (size_t j = 0; j < points[0].size(); j++) {
                glNormal3d(normals[i][j][0], normals[i][j][1], normals[i][j][2]);
                float F_schlick = schlick(specular, light_position, normals[i][j]);
                F_schlick /= 100000000000000000000.0;
                buf[0] = F_schlick;
                buf[1] = F_schlick;
                buf[2] = F_schlick;
                buf[3] = 0;
                glMaterialfv(GL_FRONT, GL_SPECULAR, buf);
                //buf[0] = 0;
                //buf[1] = 0;
                //buf[2] = 0;
                //buf[3] = 1;
                //glColor4f(F_schlick, F_schlick, F_schlick, 0.9);
                glVertex3d((i - center_x), points[i][j]->height, (j - center_z));

				glNormal3d(normals[i + 1][j][0], normals[i+1][j][1], normals[i + 1][j][2]);
                F_schlick = schlick(specular, light_position, normals[i + 1][j]);
                F_schlick /= 100000000000000000000.0;
				buf[0] = F_schlick;
				buf[1] = F_schlick;
				buf[2] = F_schlick;
                buf[3] = 0;
                glMaterialfv(GL_FRONT, GL_SPECULAR, buf);
                //glColor4f(F_schlick, F_schlick, F_schlick, 0.9);
                glVertex3d((i + 1 - center_x), points[i + 1][j]->height, (j - center_z));
            }
            glEnd();
		}
        //glDisable(GL_COLOR_MATERIAL);
    }
};

void display();
void key(unsigned char key, int x, int y);
void update(int timer_id);
void drawWaterSurface(int timer_id);
void rain(int timer_id);
void changeLightPosition(int timer_id);

WaterSurface* waterSurface;

int main(int argc, char** argv)
{
    glutInit(&argc, argv);

    waterSurface = new WaterSurface(100, 100, 0.01);
    waterSurface->material = Material(Params(0.1f, 0.1f, 0.1f, 0.1f), Params(black).alpha(1.0f), Params(black), 128);

    glutInitWindowPosition(0, 0);
    glutInitWindowSize(800, 800);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_ALPHA);
    glutCreateWindow("ìd20-0070 è¨ñxê≥é˜");
    glClearColor(0.0, 0.0, 0.0, 1.0);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, 1.0, 1.0, 50.0);
//    gluLookAt(0, 1, 0, 0, 0.5, -8, 0, 1, 0);
    gluLookAt(0, 3, 0, 0, 0, -8, 0, 1, 0);

    glShadeModel(GL_SMOOTH);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);

    // light_positionê≥ãKâª
    float len = sqrt(light_position[0] * light_position[0] + light_position[1] * light_position[1] + light_position[2] * light_position[2]);
    light_position[0] /= len;
    light_position[1] /= len;
    light_position[2] /= len;

    float pos[] = {light_position[0], light_position[1], light_position[2], 0.0};

    glLightfv(GL_LIGHT0, GL_POSITION, pos);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
    glLightfv(GL_LIGHT0, GL_SPECULAR, white);
    glLightfv(GL_LIGHT0, GL_AMBIENT, white);

    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutTimerFunc(100, update, 0);
    glutTimerFunc(100, drawWaterSurface, 1);
    glutTimerFunc(100, rain, 2);
    glutTimerFunc(100, changeLightPosition, 3);
    glutMainLoop();
	return 0;
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

    glTranslated(0, 0, -5);

    // îwåi
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
        glScaled(0.1, 0.5, 0.1);
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
    glutTimerFunc(1000 / FRAMES_PER_SECOND, drawWaterSurface, 1);
}

void update(int timer_id)
{
    if (timer_id != 0)
    {
		return;
	}
    if (waterSurface != nullptr)
    {
		waterSurface->update(0.1);
	}

	glutTimerFunc(1, update, 0);
}

void rain(int timer_id)
{
    if (timer_id != 2)
    {
		return;
	}
    
    uint64_t x = rand() % 98;
    uint64_t y = rand() % 98;

    waterSurface->addForce(x + 1, y + 1, -0.5);
    waterSurface->addForce(x + 1, y, -0.25);
    waterSurface->addForce(x + 1, y + 2, -0.25);
    waterSurface->addForce(x, y + 1, -0.25);
    waterSurface->addForce(x + 2, y + 1, -0.25);

	glutTimerFunc(50 + rand() % 100, rain, 2);
}

bool change = false;

void key(unsigned char key, int x, int y)
{
    if (key == ' ')
    {
		change = !change;
        if (change)
        {
            glutTimerFunc(50, changeLightPosition, 3);
		}
	}
}

double count = 0;

void changeLightPosition(int timer_id)
{
    if (timer_id != 3)
    {
		return;
	}

    count += 0.1;
    if (count > 2 * 3.141592)
    {
		count = 0;
    }

	float light_position[] = { 0, 1.0, 0, 0.0 };
//	light_position[0] = sin(count);
	light_position[1] = abs(cos(count));
	light_position[2] = sin(count);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

    if(change)
	glutTimerFunc(50, changeLightPosition, 3);
}