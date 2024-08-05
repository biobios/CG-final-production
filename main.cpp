#include "glut.h"
#include <GL/gl.h>
#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <Windows.h>
#include <mmsystem.h>

constexpr float FRAMES_PER_SECOND = 20;

constexpr float white[] = {1.0, 1.0, 1.0, 1.0};
constexpr float black[] = {0.0, 0.0, 0.0, 1.0};
constexpr float sunset[] = {0.92549, 0.3647059, 0.058823, 1.0};
constexpr float transparent[] = {0.0, 0.0, 0.0, 0.0};
constexpr double PI = 3.1415926535;

std::mt19937 *random_engine = nullptr;
std::uniform_real_distribution<double> random_dist(0, 1);

class Environment
{
    std::array<double, 3> eye;
    std::array<double, 3> center;
    std::array<double, 3> up;
    
    std::array<double, 3> H;

    std::array<double, 3> light_position;
    std::array<float, 4> light_color;
    void calcH();
public:
    Environment();
    void lookAt(std::array<double, 3> const& eye, std::array<double, 3> const& center, std::array<double, 3> const& up);
    void setLightPosition(std::array<double, 3> const& light_position);
    void setLightColor(std::array<float, 4> const& light_color);
    void calcSpecular(std::array<float, 4> const& material, std::array<double, 3> const& normal, double smoothness, std::array<float, 4>& output_specular);
    void apply();
    float calcFresnelReflectance(double c_spec, std::array<double, 3> const &normal);
};

class Transition
{
    std::array<double, 3> next_light_position;
    std::array<float, 4> next_light_color;

    std::array<double, 3> prev_light_position;
    std::array<float, 4> prev_light_color;

    uint64_t transition_time_ms;
    uint64_t current_time = 0;
    double prev_raindrops_per_second;
    double next_raindrops_per_second;

public:
    Transition(std::array<double, 3> const &prev_light_position, std::array<float, 4> const &prev_light_color, double prev_raindrops_per_second) : prev_light_position(prev_light_position), prev_light_color(prev_light_color), prev_raindrops_per_second(prev_raindrops_per_second)
    {
        // 0~0.7なら
        // next_light_colorをwhiteからsunsetの間でランダムに決定する
        // 0.7~1.0なら
        // next_light_colorをwhiteの0.7~1.0倍の間でランダムに決定する

        double random_value = random_dist(*random_engine);
        if (random_value < 0.7)
        {
            random_value = random_dist(*random_engine);
            this->next_light_color[0] = white[0] * (1 - random_value) + sunset[0] * random_value;
            this->next_light_color[1] = white[1] * (1 - random_value) + sunset[1] * random_value;
            this->next_light_color[2] = white[2] * (1 - random_value) + sunset[2] * random_value;
            this->next_light_color[3] = 1.0;
        }
        else
        {
            this->next_light_color[0] = white[0] * (1 - random_value) + white[0] * random_value;
            this->next_light_color[1] = white[1] * (1 - random_value) + white[1] * random_value;
            this->next_light_color[2] = white[2] * (1 - random_value) + white[2] * random_value;
            this->next_light_color[3] = 1.0;
        }

        // radをπ/8から2π/6の間でランダムに決定する
        random_value = random_dist(*random_engine);
        double rad = PI / 8 * (1 - random_value) + 2 * PI / 6 * random_value;

        this->next_light_position[1] = sin(rad);
        this->next_light_position[2] = -cos(rad);
        // radを-π/12からπ/12の間でランダムに決定する
        random_value = random_dist(*random_engine);
        rad = -PI / 12 * (1 - random_value) + PI / 12 * random_value;
        this->next_light_position[0] = sin(rad);

        // 1分から5分の間でランダムに決定する
        random_value = random_dist(*random_engine);
        //        this->transition_time_ms = 60000 * (1 - random_value) + 300000 * random_value;
        this->transition_time_ms = 15000;

        // raindrops_per_secondを2から10の間でランダムに決定する
        random_value = random_dist(*random_engine);
        this->next_raindrops_per_second = 2 * (1 - random_value) + 10 * random_value;
    }

    void applyEnvironment(Environment *env)
    {
        std::array<double, 3> current_light_position;
        std::array<float, 4> current_light_color;

        double transition_ratio = (double)current_time / transition_time_ms;

        current_light_position[0] = prev_light_position[0] * (1 - transition_ratio) + next_light_position[0] * transition_ratio;
        current_light_position[1] = prev_light_position[1] * (1 - transition_ratio) + next_light_position[1] * transition_ratio;
        current_light_position[2] = prev_light_position[2] * (1 - transition_ratio) + next_light_position[2] * transition_ratio;

        current_light_color[0] = prev_light_color[0] * (1 - transition_ratio) + next_light_color[0] * transition_ratio;
        current_light_color[1] = prev_light_color[1] * (1 - transition_ratio) + next_light_color[1] * transition_ratio;
        current_light_color[2] = prev_light_color[2] * (1 - transition_ratio) + next_light_color[2] * transition_ratio;
        current_light_color[3] = 1.0;

        env->setLightPosition(current_light_position);
        env->setLightColor(current_light_color);
    }

    void update(uint64_t dt)
    {
        current_time += dt;
    }

    uint64_t getTimeUntilNextRaindrop() const
    {
        double current_raindrops_per_second = prev_raindrops_per_second * (1 - (double)current_time / transition_time_ms) + next_raindrops_per_second * (double)current_time / transition_time_ms;
        std::exponential_distribution<double> dist(current_raindrops_per_second);
        return dist(*random_engine) * 1000;
    }

    bool isFinished() const
    {
        return current_time >= transition_time_ms;
    }

    Transition *createNextTransition()
    {
        return new Transition(next_light_position, next_light_color, next_raindrops_per_second);
    }
};

Environment *env;
Transition *transition = nullptr;

class MassPoint
{
public:
    double height;
    double velocity;
    static constexpr double RESISTANCE = 0.1;
    static constexpr double GRAVITY = 0.01;
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
        velocity += dv - velocity * RESISTANCE * dt;
        velocity -= height * GRAVITY * dt;
//        velocity += dv;
        height += velocity * dt;
    }
};

struct Params {
    std::array<float, 4> p;
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
    operator float* () { return p.data(); }
    operator std::array<float, 4>& () { return p; }
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

inline float scaled_sigmoide(double x, double gain, double mid) {
    float min = 1/ (1 + exp(-gain * (-mid)));
    float max = 1 / (1 + exp(-gain * (1 - mid)));

    return (1 / (1 + exp(-gain * (x - mid))) - min) / (max - min);
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

    void draw(float gain = 1.0f)
    {

        std::vector<std::vector<std::array<double, 3>>> normals(points.size(), std::vector<std::array<double, 3>>(points[0].size(), {gain, NORMAL_Y, gain}));
        
        for (size_t i = 1; i < points.size() - 1; i++) {
            for (size_t j = 1; j < points[0].size() - 1; j++) {
                normals[i][j][0] *= points[i + 1][j]->height - points[i - 1][j]->height;
                normals[i][j][2] *= points[i][j - 1]->height - points[i][j + 1]->height;
            }
        }

        size_t last_x = points.size() - 1;
        size_t last_y = points[0].size() - 1;

        for (size_t i = 1; i < points.size() - 1; i++) {
            normals[i][0][0] *= points[i + 1][0]->height - points[i - 1][0]->height;
            normals[i][0][2] *= (points[i][0]->height - points[i][1]->height) * 2;

            normals[i][last_y][0] *= points[i + 1][last_y]->height - points[i - 1][last_y]->height;
            normals[i][last_y][2] *= (points[i][last_y - 1]->height - points[i][last_y]->height) * 2;
        }

        for (size_t i = 1; i < points[0].size() - 1; i++) {
			normals[0][i][0] *= (points[1][i]->height - points[0][i]->height) * 2;
			normals[0][i][2] *= points[0][i - 1]->height - points[0][i + 1]->height;

			normals[last_x][i][0] *= (points[last_x][i]->height - points[last_x - 1][i]->height) * 2;
			normals[last_x][i][2] *= points[last_x][i - 1]->height - points[last_x][i + 1]->height;
        }

        normals[0][0][0] *= (points[1][0]->height - points[0][0]->height) * 2;
        normals[0][0][2] *= (points[0][0]->height - points[0][1]->height) * 2;

        normals[0][last_y][0] *= (points[1][last_y]->height - points[0][last_y]->height) * 2;
        normals[0][last_y][2] *= (points[0][last_y - 1]->height - points[0][last_y]->height) * 2;

        normals[last_x][0][0] *= (points[last_x][0]->height - points[last_x - 1][0]->height) * 2;
        normals[last_x][0][2] *= (points[last_x][0]->height - points[last_x][1]->height) * 2;

        normals[last_x][last_y][0] *= (points[last_x][last_y]->height - points[last_x - 1][last_y]->height) * 2;
        normals[last_x][last_y][2] *= (points[last_x][last_y - 1]->height - points[last_x][last_y]->height) * 2;

        // ���K��

        for (size_t i = 0; i < points.size(); i++) {
            for (size_t j = 0; j < points[0].size(); j++) {
                double len = sqrt(normals[i][j][0] * normals[i][j][0] + NORMAL_Y * NORMAL_Y + normals[i][j][2] * normals[i][j][2]);
				normals[i][j][0] /= len;
                normals[i][j][1] = NORMAL_Y / len;
				normals[i][j][2] /= len;
			}
		}
        
        // スペキュラーの計算
        std::vector<std::vector<float>> speculars(points.size(), std::vector<float>(points[0].size(), 0));

        for (size_t i = 0; i < points.size(); i++)
        {
            for (size_t j = 0; j < points[0].size(); j++)
            {
                speculars[i][j] = env->calcFresnelReflectance(0.2, normals[i][j]);
            }
        }

        double center_x = (points.size() - 1) / 2.0;
        double center_z = (points[0].size() - 1) / 2.0;

        material.apply();

        glEnable(GL_COLOR_MATERIAL);
        glEnable(GL_NORMALIZE);
        glColorMaterial(GL_FRONT, GL_SPECULAR);

        float F;

        for (size_t j = 0; j < points[0].size() - 1; j++)
        {
            glBegin(GL_TRIANGLE_STRIP);
            for (size_t i = 0; i < points.size(); i++)
            {
                glNormal3d(normals[i][j][0], normals[i][j][1], normals[i][j][2]);
                F = speculars[i][j];
                glColor4f(F, F, F, 1);
                glVertex3d((i - center_x), points[i][j]->height * gain, (j - center_z));

				glNormal3d(normals[i][j + 1][0], normals[i][j + 1][1], normals[i][j + 1][2]);
                F = speculars[i][j + 1];
                glColor4f(F, F, F, 1);
                glVertex3d((i - center_x), points[i][j + 1]->height * gain, (j + 1 - center_z));
            }
            glEnd();
		}

        glDisable(GL_COLOR_MATERIAL);
        glDisable(GL_NORMALIZE);
    }
};

// Environment(private)

void Environment::calcH()
{
    double ray[3] = { center[0] - eye[0], center[1] - eye[1], center[2] - eye[2] };
    double len = sqrt(ray[0] * ray[0] + ray[1] * ray[1] + ray[2] * ray[2]);
    ray[0] /= len;
    ray[1] /= len;
    ray[2] /= len;

    H = { light_position[0] - ray[0], light_position[1] - ray[1], light_position[2] - ray[2] };
	len = sqrt(H[0] * H[0] + H[1] * H[1] + H[2] * H[2]);
	H[0] /= len;
	H[1] /= len;
	H[2] /= len;
}

// Environment

Environment::Environment()
{
	eye = { 0, 0, 0 };
	center = { 0, 0, -1 };
	up = { 0, 1, 0 };
	light_position = { 0, 1, 0 };
	light_color = { 1, 1, 1, 1 };
	calcH();
}

void Environment::lookAt(std::array<double, 3> const& eye, std::array<double, 3> const& center, std::array<double, 3> const& up)
{
	this->eye = eye;
	this->center = center;
	this->up = up;
    
    calcH();
}

void Environment::setLightPosition(std::array<double, 3> const& light_position)
{
	this->light_position = light_position;

    double len = sqrt(light_position[0] * light_position[0] + light_position[1] * light_position[1] + light_position[2] * light_position[2]);
    this->light_position[0] /= len;
    this->light_position[1] /= len;
    this->light_position[2] /= len;

    calcH();
}

void Environment::setLightColor(std::array<float, 4> const& light_color)
{
	this->light_color = light_color;
}

void Environment::calcSpecular(std::array<float, 4> const& material, std::array<double, 3> const& normal, double smoothness, std::array<float, 4>& output_specular)
{
	double c_spec = 0.2;
	double c = light_position[0] * normal[0] + light_position[1] * normal[1] + light_position[2] * normal[2];

    double F_schlick;
    if(c >= 0)
        F_schlick = c_spec + (1 - c_spec) * pow(1 - c, 5);
    else
        F_schlick = c_spec + (1 - c_spec) * pow(1 + c, 5);
    double angle_H = acos(H[0] * normal[0] + H[1] * normal[1] + H[2] * normal[2]);
    angle_H /= smoothness;
    double k_spec = exp(- angle_H * angle_H);

	output_specular[0] = material[0] * light_color[0] * F_schlick * k_spec;
    output_specular[1] = material[1] * light_color[1] * F_schlick * k_spec;
    output_specular[2] = material[2] * light_color[2] * F_schlick * k_spec;
    output_specular[3] = material[3] * light_color[3] * F_schlick * k_spec;
//    output_specular[3] = scaled_sigmoide(output_specular[3], 100, 0.55);
}

float Environment::calcFresnelReflectance(double c_spec, std::array<double, 3> const& normal)
{
	double c = light_position[0] * normal[0] + light_position[1] * normal[1] + light_position[2] * normal[2];

    if(c >= 0)
		return c_spec + (1 - c_spec) * pow(1 - c, 5);
	else
	    return c_spec + (1 - c_spec) * pow(1 + c, 5);
}

void Environment::apply()
{


    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eye[0], eye[1], eye[2], center[0], center[1], center[2], up[0], up[1], up[2]);
    float pos[] = { light_position[0], light_position[1], light_position[2], 0.0 };
    glLightfv(GL_LIGHT0, GL_POSITION, pos);

    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_color.data());
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_color.data());
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_color.data());
}


void display();
void key(unsigned char key, int x, int y);
void update(int timer_id);
void drawWaterSurface(int timer_id);
void rain(int timer_id);
void changeLightPosition(int timer_id);
void changeEyePosition(int timer_id);
void reshapeWindow(int w, int h);

WaterSurface* waterSurface;

int main(int argc, char** argv)
{
    glutInit(&argc, argv);

    waterSurface = new WaterSurface(100, 100, 0.06);
    waterSurface->material = Material(Params(black), Params(white).maltiply(0.01f), Params(white), 128);

    env = new Environment();

    // 乱数生成器の初期化
    std::random_device seed_gen;
    random_engine = new std::mt19937(seed_gen());

    // Transitionの初期化
    transition = new Transition({0, 1, -6}, {1, 1, 1, 1}, 5);

    const int init_width = 1200;
    const int init_height = 800;

    glutInitWindowPosition(0, 0);
    glutInitWindowSize(init_width, init_height);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_ALPHA);
    glutCreateWindow("電20-0070 小堀正樹");
    glClearColor(0.0, 0.0, 0.0, 1.0);

    //glMatrixMode(GL_PROJECTION);
    //glLoadIdentity();
    //gluPerspective(45, 1.0, 1.0, 50.0);
//    gluLookAt(0, 1, 0, 0, 0.5, -8, 0, 1, 0);
    //gluLookAt(0, 3, 0, 0, 0, -8, 0, 1, 0);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, 1.0, 1.0, 50.0);

    env->lookAt({ 0, 5, 8 }, { 0, 0, 0 }, { 0, 1, 0 });

    glShadeModel(GL_SMOOTH);
    glBlendFunc(GL_ONE, GL_ONE);
    glEnable(GL_BLEND);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);

    env->setLightPosition({0, 1, -6});
    env->setLightColor({1, 0.5, 0.5, 1});

    // 雨音ループ再生
    PlaySound(TEXT("rain.wav"), NULL, SND_FILENAME | SND_ASYNC | SND_LOOP);

    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutTimerFunc(100, update, 0);
    glutTimerFunc(100, drawWaterSurface, 1);
    glutTimerFunc(100, rain, 2);
    glutMainLoop();
	return 0;
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    env->apply();

	glMatrixMode(GL_MODELVIEW);

    // �w�i
    glMaterialfv(GL_FRONT, GL_DIFFUSE, Params(white).maltiply(0.5));
    glMaterialfv(GL_FRONT, GL_SPECULAR, black);
    glMaterialfv(GL_FRONT, GL_AMBIENT, Params(white).maltiply(0.3));
    glMaterialf(GL_FRONT, GL_SHININESS, 0);

    const double outer_r = 30;
    const double inner_r = 10;

    // 底
    glBegin(GL_QUADS);
    glNormal3d(0, 1, 0);
    glVertex3d(-inner_r, -1, -inner_r);
    glVertex3d(-inner_r, -1, inner_r);
    glVertex3d(inner_r, -1, inner_r);
    glVertex3d(inner_r, -1, -inner_r);
    glEnd();
    // 奥壁
    glBegin(GL_QUADS);
    glNormal3d(0, 0, 1);
    glVertex3d(-inner_r, -1, -inner_r);
    glVertex3d(inner_r, -1, -inner_r);
    glVertex3d(inner_r, 0.2, -inner_r);
    glVertex3d(-inner_r, 0.2, -inner_r);
    glEnd();
    // 左壁
    glBegin(GL_QUADS);
    glNormal3d(1, 0, 0);
    glVertex3d(-inner_r, -1, -inner_r);
    glVertex3d(-inner_r, -1, inner_r);
    glVertex3d(-inner_r, 0.2, inner_r);
    glVertex3d(-inner_r, 0.2, -inner_r);
    glEnd();
    // 右壁
    glBegin(GL_QUADS);
    glNormal3d(-1, 0, 0);
    glVertex3d(inner_r, -1, -inner_r);
    glVertex3d(inner_r, -1, inner_r);
    glVertex3d(inner_r, 0.2, inner_r);
    glVertex3d(inner_r, 0.2, -inner_r);
    glEnd();
    // 奥床
    glBegin(GL_QUADS);
    glNormal3d(0, 1, 0);
    glVertex3d(-outer_r, 0.2, -inner_r);
    glVertex3d(-outer_r, 0.2, -outer_r);
    glVertex3d(outer_r, 0.2, -outer_r);
    glVertex3d(outer_r, 0.2, -inner_r);
    glEnd();
    // 左床
    glBegin(GL_QUADS);
    glNormal3d(0, 1, 0);
    glVertex3d(-outer_r, 0.2, -inner_r);
    glVertex3d(-outer_r, 0.2, inner_r);
    glVertex3d(-inner_r, 0.2, inner_r);
    glVertex3d(-inner_r, 0.2, -inner_r);
    glEnd();
    // 右床
    glBegin(GL_QUADS);
    glNormal3d(0, 1, 0);
    glVertex3d(outer_r, 0.2, -inner_r);
    glVertex3d(outer_r, 0.2, inner_r);
    glVertex3d(inner_r, 0.2, inner_r);
    glVertex3d(inner_r, 0.2, -inner_r);
    glEnd();

    if (waterSurface != nullptr)
    {
        glPushMatrix();
        glScaled(0.1, 0.1, 0.1);
        glTranslated(0, 1, 0);
        waterSurface->draw(5.0f);
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

    transition->update(1000 / FRAMES_PER_SECOND);
    transition->applyEnvironment(env);

    if (transition->isFinished())
    {
        Transition *next_transition = transition->createNextTransition();
        delete transition;
        transition = next_transition;
    }

    glutPostRedisplay();
    glutTimerFunc(1000 / FRAMES_PER_SECOND, drawWaterSurface, 1);
}

void reshapeWindow(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(30, (double)w / h, 1.0, 50.0);
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

    const int64_t r = 2;
    
    uint64_t x = rand() % (100 - r * 2) + r;
    uint64_t y = rand() % (100 - r * 2) + r;

    //waterSurface->addForce(x, y, -1);
    //waterSurface->addForce(x + 2, y, 0.1);
    //waterSurface->addForce(x - 2, y, 0.1);
    //waterSurface->addForce(x, y + 2, 0.1);
    //waterSurface->addForce(x, y - 2, 0.1);
    //waterSurface->addForce(x + 1, y + 1, 0.1);
    //waterSurface->addForce(x - 1, y - 1, 0.1);
    //waterSurface->addForce(x + 1, y - 1, 0.1);
    //waterSurface->addForce(x - 1, y + 1, 0.1);

    for (int i = -r; i <= r; i++)
    {
        for (int j = -r; j <= r; j++)
        {
            if (i * i + j * j <= r * r)
            {
				waterSurface->addForce(x + i, y + j, - 1 * exp(- (i * i + j * j) / 4.0));
                //waterSurface->addForce(x + i, y + j, 0.1 / r * (i * i + j * j) - 0.05 / r);
			}
		}
	}

    //waterSurface->addForce(rand() % 100, rand() % 100, -2);

    glutTimerFunc(transition->getTimeUntilNextRaindrop(), rain, 2);
}

bool change = false;
bool changeEye = false;

double angle_y = 0;
double angle_xz = 3.1415926535 / 4;

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
 //   if (key == 'e')
 //   {
	//	changeEye = !changeEye;
 //       if (changeEye)
 //       {
	//		glutTimerFunc(50, changeEyePosition, 4);
	//	}
	//}


    bool change_angle = false;

    if (key == 'w')
    {
        // ���_����Ɉړ�
        angle_xz += 0.1;
	    change_angle = true;
    }
    if (key == 's')
    {
		// ���_�����Ɉړ�
		angle_xz -= 0.1;
        change_angle = true;
    }
    if (key == 'a')
    {
        // ���_�����Ɉړ�
        angle_y -= 0.1;
        change_angle = true;
    }
    if (key == 'd')
    {
		// ���_���E�Ɉړ�
		angle_y += 0.1;
        change_angle = true;
	}

    // angle�𐧌�
    if(angle_xz > 3.141592 / 2)
	{
        angle_xz = 3.141592 / 2;
    }
    if (angle_xz < -3.141592 / 2)
    {
		angle_xz = -3.141592 / 2;
	}

    // angle_y��0����2�΂͈̔͂Ɏ��߂�
    if (angle_y > 2 * 3.141592)
    {
        angle_y -= 2 * 3.141592;
    }
    if (angle_y < 0)
    {
        angle_y += 2 * 3.141592;
    }
    
    if (change_angle) {
        std::array<double, 3> eye = {};
        eye[0] = sin(angle_y) * cos(angle_xz);
        eye[1] = abs(sin(angle_xz));
        eye[2] = cos(angle_y) * cos(angle_xz);

        eye[0] *= 9;
        eye[1] *= 9;
        eye[2] *= 9;

        env->lookAt(eye, { 0, 0, 0 }, { 0, 1, 0 });
    }
}

double count = 0;

void changeLightPosition(int timer_id)
{
    if (timer_id != 3)
    {
		return;
	}

    count += 0.01;
    if (count > 2 * 3.141592)
    {
		count = 0;
    }

	float light_position[] = { 0, 1.0, 0, 0.0 };
//	light_position[0] = sin(count);
	light_position[1] = abs(cos(count));
	light_position[2] = sin(count);

	env->setLightPosition({ light_position[0], light_position[1], light_position[2] });

    if(change)
	glutTimerFunc(50, changeLightPosition, 3);
    else
        std::cout << "y = " << light_position[1] << " z = " << light_position[2] << std::endl;
}

double angle = 0;

void changeEyePosition(int timer_id)
{
    if (timer_id != 4)
    {
		return;
	}

    angle += 0.005;
    if (angle > 2 * 3.141592)
    {
        angle = 0;
	}

    float eye_position[] = { 0, 1.0, 0, 0.0 };
	eye_position[0] = 3 * sin(angle);
	eye_position[1] = 9 * abs(cos(angle));
	eye_position[2] = 3 * sin(angle * 2);

	env->lookAt({ eye_position[0], eye_position[1], eye_position[2] }, { 0, 0, 0 }, { 0, 1, 0 });

	if(changeEye)
	glutTimerFunc(50, changeEyePosition, 4);
	else
		std::cout << "x = " << eye_position[0] << " y = " << eye_position[1] << " z = " << eye_position[2] << std::endl;
}