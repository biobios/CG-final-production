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

// Environment
// 視点の位置、光源の位置、光源の色を管理するクラス
class Environment
{
    std::array<double, 3> eye;
    std::array<double, 3> center;
    std::array<double, 3> up;

    std::array<double, 3> light_position;
    std::array<float, 4> light_color;

public:
    Environment();
    // 視点の位置、注視点の位置、上方向のベクトルを設定する
    void lookAt(std::array<double, 3> const &eye, std::array<double, 3> const &center, std::array<double, 3> const &up);
    // 光源の位置を設定する
    void setLightPosition(std::array<double, 3> const &light_position);
    // 光源の色を設定する
    void setLightColor(std::array<float, 4> const &light_color);
    // 環境を適用する
    void apply();
    // Fresnel反射率を計算する
    float calcFresnelReflectance(double c_spec, std::array<double, 3> const &normal);
};

// Transition
// 環境の遷移を管理するクラス
// 管理する対象には雨の落ちる頻度を含む
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
        // それをもとに, y座標とz座標を決定する
        random_value = random_dist(*random_engine);
        double rad = PI / 8 * (1 - random_value) + 2 * PI / 6 * random_value;

        this->next_light_position[1] = sin(rad);
        this->next_light_position[2] = -cos(rad);

        // radを-π/12からπ/12の間でランダムに決定する
        // それをもとに, x座標を決定する
        random_value = random_dist(*random_engine);
        rad = -PI / 12 * (1 - random_value) + PI / 12 * random_value;
        this->next_light_position[0] = sin(rad);

        // 1分から3分の間で遷移にかかる時間をランダムに決定する
        random_value = random_dist(*random_engine);
        this->transition_time_ms = 60000 * (1 - random_value) + 180000 * random_value;
        
        // raindrops_per_secondを2から10の間でランダムに決定する
        random_value = random_dist(*random_engine);
        this->next_raindrops_per_second = 2 * (1 - random_value) + 10 * random_value;
    }

    // 現在の光源の位置、光源の色を環境に適用する
    void applyEnvironment(Environment *env)
    {
        // 経過時間に応じて光源の位置、光源の色を変化させる
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

    // 経過時間を更新する
    void update(uint64_t dt)
    {
        current_time += dt;
        if (current_time > transition_time_ms)
        {
			current_time = transition_time_ms;
		}
    }

    // 次の雨粒が落ちるまでの時間を取得する
    uint64_t getTimeUntilNextRaindrop() const
    {
        // 経過時間に応じて雨粒の落ちる頻度を変化させる
        double current_raindrops_per_second = prev_raindrops_per_second * (1 - (double)current_time / transition_time_ms) + next_raindrops_per_second * (double)current_time / transition_time_ms;
        
        // 指数分布に従って次の雨粒が落ちるまでの時間をランダムに決定する
        std::exponential_distribution<double> dist(current_raindrops_per_second);
        return dist(*random_engine) * 1000;
    }

    // 遷移が終了したかどうかを判定する
    bool isFinished() const
    {
        return current_time >= transition_time_ms;
    }

    // 次の遷移を生成する
    Transition *createNextTransition()
    {
        // 現在の遷移が修了していたら、終了時の状態から次の遷移を生成する
        if(isFinished())
            return new Transition(next_light_position, next_light_color, next_raindrops_per_second);

        // 現在の遷移が修了していない場合は、現在の状態から次の遷移を生成する
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

        double current_raindrops_per_second = prev_raindrops_per_second * (1 - (double)current_time / transition_time_ms) + next_raindrops_per_second * (double)current_time / transition_time_ms;

        return new Transition(current_light_position, current_light_color, current_raindrops_per_second);

    }
};

Environment *env = nullptr;
Transition *transition = nullptr;

// MassPoint
// 質点の位置、速度を管理するクラス
class MassPoint
{
public:
    double height;
    double velocity;

    // 速度に比例する抵抗係数
    static constexpr double RESISTANCE = 0.1;
    // 原点を基準としたばね定数
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

    // 質点の位置、速度を更新する
    void update(double f, double dt)
    {
        velocity += f * dt - velocity * RESISTANCE * dt;
        velocity -= height * GRAVITY * dt;
        height += velocity * dt;
    }
};

// Params
// パラメータを管理する構造体
struct Params
{
    std::array<float, 4> p;
    Params() = default;
    Params(float r, float g, float b, float a)
    {
        p[0] = r;
        p[1] = g;
        p[2] = b;
        p[3] = a;
    }
    Params(float const *p)
    {
        for (int i = 0; i < 4; i++)
        {
            this->p[i] = p[i];
        }
    }

    // 0~2番目の要素にaを掛ける
    Params &maltiply(float a)
    {
        for (int i = 0; i < 3; i++)
        {
            p[i] *= a;
        }
        return *this;
    }

    // 3番目の要素にaを代入する
    Params &alpha(float a)
    {
        p[3] = a;
        return *this;
    }

    // ポインタへの変換関数
    operator float *() { return p.data(); }
    // arrayへの変換関数
    operator std::array<float, 4> &() { return p; }
};

// Material
// マテリアルを管理する構造体
struct Material
{
    Params ambient;
    Params diffuse;
    Params specular;
    float shininess;
    Material() = default;
    Material(Params const &a, Params const &d, Params const &s, float sh) : ambient(a), diffuse(d), specular(s), shininess(sh) {}
    // マテリアルを適用する
    void apply()
    {
        glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
        glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
        glMaterialf(GL_FRONT, GL_SHININESS, shininess);
    }
};

// WaterSurface
// 水面を管理するクラス
class WaterSurface
{
    std::vector<std::vector<MassPoint *>> points;
    double tension;
    double specular;

public:
    Material material = {};

    WaterSurface(uint64_t n, uint64_t m, double t, double s = 0.02) : tension(t), specular(s)
    {
        // n*mの質点を生成する
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

    // 質点に力を加える
    void addForce(uint64_t i, uint64_t j, double f)
    {
        points[i][j]->addVelocity(f);
    }

    // 質点の高さを設定する
    void setHeight(uint64_t i, uint64_t j, double h)
    {
        points[i][j]->height = h;
    }

    // dt秒後の状態に更新する
    void update(double dt)
    {
        // 各質点について、隣接する質点とのy座標の差の和を計算する
        std::vector<std::vector<double>> force(points.size(), std::vector<double>(points[0].size(), 0));

        // x方向の隣接する質点とのy座標の差の和を計算する
        for (size_t i = 0; i < points.size() - 1; i++)
        {
            for (size_t j = 0; j < points[0].size(); j++)
            {
                force[i][j] += points[i + 1][j]->height;
                force[i + 1][j] += points[i][j]->height;

                force[i][j] -= points[i][j]->height;
                force[i + 1][j] -= points[i + 1][j]->height;
            }
        }

        // y方向の隣接する質点とのy座標の差の和を計算する
        for (size_t i = 0; i < points.size(); i++)
        {
            for (size_t j = 0; j < points[0].size() - 1; j++)
            {
                force[i][j] += points[i][j + 1]->height;
                force[i][j + 1] += points[i][j]->height;

                force[i][j] -= points[i][j]->height;
                force[i][j + 1] -= points[i][j + 1]->height;
            }
        }

        // 質点の位置を更新する
        for (size_t i = 0; i < points.size(); i++)
        {
            for (size_t j = 0; j < points[0].size(); j++)
            {
                points[i][j]->update(force[i][j] * tension, dt);
            }
        }
    }

    static constexpr double NORMAL_Y = 2;

    void draw(float gain = 1.0f)
    {
        // 法線ベクトルの計算
        // 質点の位置関係を上から見て次のように表す
        //   2
        // 3 0 1
        //   4
        // この時、質点0の法線ベクトルを、
        // (1 - 3) × (2 - 4)
        // として計算すると
        // n = (y3 - y1, 2, y2 - y4)となる
        std::vector<std::vector<std::array<double, 3>>> normals(points.size(), std::vector<std::array<double, 3>>(points[0].size(), {gain, NORMAL_Y, gain}));

        // 端を除いた部分(全方向に隣が存在する部分)の法線ベクトルを計算する
        for (size_t i = 1; i < points.size() - 1; i++)
        {
            for (size_t j = 1; j < points[0].size() - 1; j++)
            {
                normals[i][j][0] *= points[i + 1][j]->height - points[i - 1][j]->height;
                normals[i][j][2] *= points[i][j - 1]->height - points[i][j + 1]->height;
            }
        }

        size_t last_x = points.size() - 1;
        size_t last_y = points[0].size() - 1;

        // 角を除く左右の端の法線ベクトルを計算する
        for (size_t i = 1; i < points.size() - 1; i++)
        {
            normals[i][0][0] *= points[i + 1][0]->height - points[i - 1][0]->height;
            normals[i][0][2] *= (points[i][0]->height - points[i][1]->height) * 2;

            normals[i][last_y][0] *= points[i + 1][last_y]->height - points[i - 1][last_y]->height;
            normals[i][last_y][2] *= (points[i][last_y - 1]->height - points[i][last_y]->height) * 2;
        }

        // 角を除く上下の端の法線ベクトルを計算する
        for (size_t i = 1; i < points[0].size() - 1; i++)
        {
            normals[0][i][0] *= (points[1][i]->height - points[0][i]->height) * 2;
            normals[0][i][2] *= points[0][i - 1]->height - points[0][i + 1]->height;

            normals[last_x][i][0] *= (points[last_x][i]->height - points[last_x - 1][i]->height) * 2;
            normals[last_x][i][2] *= points[last_x][i - 1]->height - points[last_x][i + 1]->height;
        }

        // 角の法線ベクトルを計算する
        normals[0][0][0] *= (points[1][0]->height - points[0][0]->height) * 2;
        normals[0][0][2] *= (points[0][0]->height - points[0][1]->height) * 2;

        normals[0][last_y][0] *= (points[1][last_y]->height - points[0][last_y]->height) * 2;
        normals[0][last_y][2] *= (points[0][last_y - 1]->height - points[0][last_y]->height) * 2;

        normals[last_x][0][0] *= (points[last_x][0]->height - points[last_x - 1][0]->height) * 2;
        normals[last_x][0][2] *= (points[last_x][0]->height - points[last_x][1]->height) * 2;

        normals[last_x][last_y][0] *= (points[last_x][last_y]->height - points[last_x - 1][last_y]->height) * 2;
        normals[last_x][last_y][2] *= (points[last_x][last_y - 1]->height - points[last_x][last_y]->height) * 2;

        // 法線ベクトルを正規化する
        for (size_t i = 0; i < points.size(); i++)
        {
            for (size_t j = 0; j < points[0].size(); j++)
            {
                double len = sqrt(normals[i][j][0] * normals[i][j][0] + NORMAL_Y * NORMAL_Y + normals[i][j][2] * normals[i][j][2]);
                normals[i][j][0] /= len;
                normals[i][j][1] = NORMAL_Y / len;
                normals[i][j][2] /= len;
            }
        }

        // 法線ベクトルからFresnel反射率を計算する
        std::vector<std::vector<float>> speculars(points.size(), std::vector<float>(points[0].size(), 0));

        for (size_t i = 0; i < points.size(); i++)
        {
            for (size_t j = 0; j < points[0].size(); j++)
            {
                speculars[i][j] = env->calcFresnelReflectance(0.2, normals[i][j]);
            }
        }

        // 中央が0,0になるように描画する
        double center_x = (points.size() - 1) / 2.0;
        double center_z = (points[0].size() - 1) / 2.0;

        material.apply();

        // glColor4fで指定した色を反射率として使用するように設定する
        glEnable(GL_COLOR_MATERIAL);
        glColorMaterial(GL_FRONT, GL_SPECULAR);

        // 法線ベクトルの正規化を有効にする(scaleの影響を受けないようにする)
        glEnable(GL_NORMALIZE);

        // 水面を描画する
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

        // 各種設定を元に戻す
        glDisable(GL_COLOR_MATERIAL);
        glDisable(GL_NORMALIZE);
    }
};

// Environment 実装

Environment::Environment()
{
    eye = {0, 0, 0};
    center = {0, 0, -1};
    up = {0, 1, 0};
    light_position = {0, 1, 0};
    light_color = {1, 1, 1, 1};
}

void Environment::lookAt(std::array<double, 3> const &eye, std::array<double, 3> const &center, std::array<double, 3> const &up)
{
    this->eye = eye;
    this->center = center;
    this->up = up;
}

void Environment::setLightPosition(std::array<double, 3> const &light_position)
{
    this->light_position = light_position;

    // 光源の位置を正規化する
    double len = sqrt(light_position[0] * light_position[0] + light_position[1] * light_position[1] + light_position[2] * light_position[2]);
    this->light_position[0] /= len;
    this->light_position[1] /= len;
    this->light_position[2] /= len;
}

void Environment::setLightColor(std::array<float, 4> const &light_color)
{
    this->light_color = light_color;
}

float Environment::calcFresnelReflectance(double c_spec, std::array<double, 3> const &normal)
{
    // 光源の位置と法線ベクトルの内積を計算する
    double c = light_position[0] * normal[0] + light_position[1] * normal[1] + light_position[2] * normal[2];

    // Fresnel反射率を計算する
    if (c >= 0)
        return c_spec + (1 - c_spec) * pow(1 - c, 5);
    else
        return c_spec + (1 - c_spec) * pow(1 + c, 5);
}

// 環境を適用する
void Environment::apply()
{
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eye[0], eye[1], eye[2], center[0], center[1], center[2], up[0], up[1], up[2]);
    float pos[] = {light_position[0], light_position[1], light_position[2], 0.0};
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
void reshapeWindow(int w, int h);

WaterSurface *waterSurface;

int main(int argc, char **argv)
{
    glutInit(&argc, argv);

    // 水面の初期化
    waterSurface = new WaterSurface(100, 100, 0.1);
    waterSurface->material = Material(Params(black), Params(white).maltiply(0.01f), Params(white), 128);

    // 環境の初期化
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

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, (double)init_width / init_height, 1.0, 50.0);

    // 視点の設定
    const double r = 16;
    std::array<double, 3> eye = {};
    eye[0] = 0;
    eye[1] = sin(PI / 6);
    eye[2] = cos(PI / 6);

    eye[0] *= r;
    eye[1] *= r;
    eye[2] *= r;
    env->lookAt(eye, {0, 0, 0}, {0, 1, 0});

    // アルファブレンドの設定
    glBlendFunc(GL_ONE, GL_ONE);
    glEnable(GL_BLEND);

    // ライティングの設定
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);

    // ライトの設定
    env->setLightPosition({0, 1, -6});
    env->setLightColor({1, 0.5, 0.5, 1});

    // 雨音ループ再生
    PlaySound(TEXT("rain.wav"), NULL, SND_FILENAME | SND_ASYNC | SND_LOOP);

    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutTimerFunc(100, update, 0);
    glutTimerFunc(100, drawWaterSurface, 1);
    glutTimerFunc(100, rain, 2);
    glutReshapeFunc(reshapeWindow);
    glutMainLoop();
    return 0;
}

// ディスプレイハンドラ
void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    env->apply();
    glMatrixMode(GL_MODELVIEW);

    glTranslated(0, 0, -4);

    // 背景の描画
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
        glScaled(0.2, 0.2, 0.2);
        waterSurface->draw(7.5f);
        glPopMatrix();
    }

    glFlush();
}

// アニメーションが行われるように一定周期で呼び出される関数
void drawWaterSurface(int timer_id)
{
    if (timer_id != 1)
    {
        return;
    }

    // 遷移に経過時間を適用する
    transition->update(1000 / FRAMES_PER_SECOND);
    transition->applyEnvironment(env);

    // 遷移が終了したら次の遷移を生成する
    if (transition->isFinished())
    {
        Transition *next_transition = transition->createNextTransition();
        delete transition;
        transition = next_transition;
    }

    // レンダリング
    glutPostRedisplay();
    glutTimerFunc(1000 / FRAMES_PER_SECOND, drawWaterSurface, 1);
}

// ウィンドウサイズが変更されたときに呼び出される関数
void reshapeWindow(int w, int h)
{
    // 表示がゆがまないように設定する
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(30, (double)w / h, 1.0, 50.0);
}

// シミュレーションが行われるように一定周期で呼び出される関数
void update(int timer_id)
{
    if (timer_id != 0)
    {
        return;
    }
    if (waterSurface != nullptr)
    {
        waterSurface->update(0.3);
    }

    glutTimerFunc(4, update, 0);
}

// 雨粒が落ちるように不定期で呼び出される関数
void rain(int timer_id)
{
    if (timer_id != 2)
    {
        return;
    }

    // 雨粒の落下サイズ 2r * 2r
    const int64_t r = 2;

    // 位置をランダムに決定
    uint64_t x = rand() % (100 - r * 2) + r;
    uint64_t y = rand() % (100 - r * 2) + r;

    for (int i = -r; i <= r; i++)
    {
        for (int j = -r; j <= r; j++)
        {
            if (i * i + j * j <= r * r)
            {
                // 中心からの距離に応じて力を加える
                waterSurface->addForce(x + i, y + j, -1 * exp(-(i * i + j * j) / 4.0));
            }
        }
    }

    // 次の雨粒が落ちるまでの時間を取得して、その時間後に再度この関数を呼び出す
    glutTimerFunc(transition->getTimeUntilNextRaindrop(), rain, 2);
}

void key(unsigned char key, int x, int y)
{
    // スペースキーを押すと遷移の終了を待たずに次の遷移に移る
    if (key == ' ')
    {
        Transition *next_transition = transition->createNextTransition();
        delete transition;
        transition = next_transition;
    }

    // qキーを押すとプログラムを終了する
    if (key == 'q')
    {
		exit(0);
	}
}