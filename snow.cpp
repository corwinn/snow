/**** BEGIN LICENSE BLOCK ****

BSD 3-Clause License

Copyright (c) 2023, the wind.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

**** END LICENCE BLOCK ****/

#ifdef _WIN32
#include <windows.h>
#endif
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include <GL/gl.h>
#include <GL/glut.h>
#include <GL/freeglut.h>
#ifdef _WIN32
#include <GL/glext.h>
#endif

#define private private:
#define public public:

template <typename T> class DMem final
{
    private T * _p;
    public DMem(size_t n)
    {
        _p = static_cast<T *>(calloc (n, sizeof(T)));
        if (! _p) printf (
            "Failed to allocate %zu bytes of memory. Stop.", n * sizeof(T)),
            exit (1);
    }
    public ~DMem() { if (_p) free (_p), _p = nullptr; }
    operator T*() { return _p; }

    public DMem(const DMem<T> &) = delete;
    public DMem operator=(const DMem<T> &) = delete;
    public DMem(DMem<T> && s) { operator= (s); }
    public DMem<T> && operator=(DMem<T> && s)
    {
        return _p = s._p, s._p = nullptr, static_cast<DMem<T> &&>(*this);
    }
};

// glFrustum parameters; _WIN32 doesn't like NEAR and FAR
static float const NEARR = 1.f, FARR = 3.f;//, DEPTH = FAR - NEAR;
static float const LEFT = -1.f, RIGHT = 1.f, BOTTOM = -1.f, TOP = 1.f;
static bool paused = false;
static GLsizei current_viewport_width {0}, current_viewport_height {0};

// float should be enough
#define float float
#define MY_PI 3.14159265358979323846
//TODO os.h
#ifdef _WIN32
float rndf(float n) { return rand () / static_cast<float>(RAND_MAX) * n; }
#else
float rndf(float n) { return random () / static_cast<float>(RAND_MAX) * n; }
#endif
float rndfr(float a, float b) { return a + rndf (b-a); }

static float const TMAX {1.0/32}; // translate max value
static float const TMIN {TMAX/10}; // translate min value
static float const RMAX {2.0}; // rotate max value
static float const RMIN {-2.0}; // rotate min value

class Snowflake final
{
    static int const TNUM {4};
    GLfloat _color;  // R=G=B
    float _wx, _wy, _wz;     // weight (0;1] by which _d are multiplied
    float _p[3*(TNUM+1)] {}; // base pos x,y,z + trail (TNUM * x,y,z)
    float _r[3*(TNUM+1)] {}; // base rotation x,y,z + trail (TNUM * x,y,z)
    int   _mblur_cnt {};     // how many of the "trail" are initialzied
    float _dx, _dy, _dz;     // delta x, y, z - translate
    float _rdx, _rdy, _rdz;  // delta x, y, z - rotate
                             // when -1 or +1 is reached, invert the _rd
    // float _ra {};            // rotation angle
    float _scale;
    float _color_fade_out {}; // _color delta for fading out; 2^n
    bool _fade_out {};
    const GLfloat _v[3*4]; // vertices: 0-3 = tl, bl, br, tr
    const GLfloat _t[2*4]; // uv
    const GLubyte _i[4] {0, 1, 3, 2}; // indices; TRISTRIP
#ifdef _WIN32
    struct rnd_init { rnd_init() { srand (time (nullptr)); } };
#else
    struct rnd_init { rnd_init() { srandom (time (nullptr)); } };
#endif
    // bounding box
    private float
        min_x{-0.25}, max_x{1.25}, min_y{-0.5}, max_y{1.5}, min_z{}, max_z{1};
    private void Init()
    {
        // cleanup the prev. mblur
        memset (_p+3, 0, 3*TNUM*sizeof(float));
        memset (_r, 0, 3*(TNUM+1)*sizeof(float));
        _mblur_cnt = 0;
        _fade_out = false;
        _wx = _wy = _wz = 1;
        _dx = rndfr (TMIN, TMAX) - rndfr (TMIN, TMAX);
        _dy = -rndfr (TMIN, TMAX);
        _dz = (rndfr (TMIN, TMAX) - rndfr (TMIN, TMAX))/10.0;
        _rdx = rndfr (RMIN, RMAX) - rndfr (RMIN, RMAX);
        _rdy = rndfr (RMIN, RMAX) - rndfr (RMIN, RMAX);
        _rdz = rndfr (RMIN, RMAX) - rndfr (RMIN, RMAX);
        _scale = rndfr (TMIN, TMAX);
        _color = rndfr (0.5, 1.0);
        _color_fade_out = {};
        _p[0] = rndfr (min_x, max_x); _p[1] = rndf (max_y); _p[2] = rndf (max_z);
        // printf ("Snowflake: pos(%5.2f,%5.2f,%5.2f)\n", _p[0], _p[1], _p[2]);
        // printf ("Snowflake: ro(%5.2f,%5.2f,%5.2f)\n", _rdx, _rdy, _rdz);
    }
    public Snowflake()
        : _v {0,1,0, 0,0,0, 1,0,0, 1,1,0}, _t {0,1, 0,0, 1,0, 1,1}
    {
        static rnd_init __r__i__ {};
        Init ();
    }
    public void Render()//TODO shader (as an option)
    {
        if (_fade_out) {
            _color -= _color_fade_out;
            if (_color < 0) return;
            _color_fade_out += _color_fade_out;
        }
        for (int t = _mblur_cnt; t >= 0; t--) {
            glLoadIdentity ();

            glTranslatef ((_p+3*t)[0], (_p+3*t)[1], (_p+3*t)[2]);
            glTranslatef (0.5, 0.5, 0);
            // glRotatef (_ra, _rx, _ry, _rz);

            glRotatef ((_r+3*t)[0], 1, 0, 0);
            glRotatef ((_r+3*t)[1], 0, 1, 0);
            glRotatef ((_r+3*t)[2], 0, 0, 1);

            glScalef (_scale, _scale, 1);
            glTranslatef (-0.5, -0.5, 0);

            glColor4f (_color, _color, _color, 1.0f / (t*t*t+1));
            glBegin (GL_TRIANGLE_STRIP);//TODO come on; (leave this optional)
                for (int i = 0; i < 4; i++)
                    glTexCoord2fv (&(_t[_i[i]*2])),
                    glVertex3fv (&(_v[_i[i]*3]));
            glEnd ();
        }
    }
    private inline Snowflake & Rotate(float & a, float & d)
    {
        a += d; if (a > 360) a = a-360; else if (a < -360) a = -a-360;
        return * this;
    }
    public void Step()
    {
        if (_color < 0) Init ();
        memmove (_p+3, _p, 3*TNUM*sizeof(float));
        memmove (_r+3, _r, 3*TNUM*sizeof(float));
        if (_mblur_cnt < TNUM) _mblur_cnt++;
        Rotate (_r[0], _rdx).Rotate (_r[1], _rdy).Rotate (_r[2], _rdz);
        _p[0] += _wx * _dx;
        _p[1] += _wy * _dy;
        _p[2] += _wz * _dz;
        // stop disappearing out of thin air
        if (_p[2] > max_z || _p[2] < min_z) _dz = -_dz, _p[2] += _dz;
        if (_fade_out) return;
        _fade_out = _p[1] < min_y
            || _p[0] < min_x || _p[0] > max_x;
        if (_fade_out) _color_fade_out = _color / 256; // ~8 frames
    }
};// Snowflake

static int const WT_SIZE {100};
static int const WT_NUM {WT_SIZE*WT_SIZE};

static float WTable[WT_NUM] {};// Shall be used by the wx, wy, wz above

static float Noise(int i, int s, float & m, float & n)
{
    // this is the shortest (LOC) ring buffer I've ever created
    int const R = 1*WT_SIZE; // MA(R)
    static float ring[R] {}, sum {};
    static int ring_p {R-1}, ring_g {0}; // put and get ptr
    // static float t {}, t2 {0.01}, t2_p {1};
#ifdef _WIN32
    float r = rand () / static_cast<float>(RAND_MAX) * s;
#else
    float r = random () / static_cast<float>(RAND_MAX) * s;
#endif
    sum -= ring[ring_g];
    sum += ring[ring_p] = r;
    ring_p = (ring_p + 1) % R;
    ring_g = (ring_g + 1) % R;
    r = sum / R;

    int x = i % WT_SIZE, y = i / WT_SIZE;
    float tx = (x+1) / static_cast<float>(WT_SIZE),
        ty = (y+1) / static_cast<float>(WT_SIZE);

    //TODO variate these coefficients to create animation
    tx *= 2.2/r,
        ty *= 23.2/r; // shift randomly

    // looks great without the "if (tx > 0.5)" above
    float y_f = WT_SIZE / 2.0, x_f = y_f / 2.0;
    float rr = tx*cos(x/x_f)*r - ty*cos(y/y_f)*r; // looks great

    if (rr > m) m = rr;
    if (rr < n) n = rr;
    return rr;
}

static class Snow final
{
    static int const NUM {200};
    Snowflake _s[NUM];
    public void Render() { for (auto & s : _s) s.Render (); }
    public void Step() { for (auto & s : _s) s.Step (); }
} Snow;

static void glutDisplay()
{
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    Snow.Render ();

    //TODO this axis isn't scaled by the projection matrix - I'm missing
    //     something obvious it seems?
    /*glLoadIdentity ();
    // glRotatef (10,1,0,0); glRotatef (-10,0,1,0);
    GLfloat a[6] {0,0,0,1,0,0};
    glBegin (GL_LINES);
        for (int i = 1; i < 4; i++)
            glColor3fv (a+i), glVertex3fv (a+0), glVertex3fv (a+i);
        glColor3fv (a+2); glVertex3fv (a+1), glVertex3f (0,0.1,1);
        glColor3fv (a+3); glVertex3fv (a+1), glVertex3f (0.1,0,1);
    glEnd ();*/

    glutSwapBuffers ();
}

// auto-hide the mouse cursor on inactivity when full-screen
static time_t mouse_moved, tmouse;
static bool fs = false;

static void glutIdle()
{
    if (fs) {
        auto df = difftime (tmouse = time (0), mouse_moved);
        if (df > 3.0)
            glutSetCursor (GLUT_CURSOR_NONE);
    }
    if (! paused) Snow.Step ();
    usleep (40000), glutPostRedisplay (); //TODO timing
}

static void glutReshape(GLsizei width, GLsizei height)
{
    if (height < 1) height = 1;
    if (width < 1) width = 1;
    current_viewport_height = height;
    current_viewport_width = width;
    glViewport (0, 0, width, height);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    glFrustum (LEFT, RIGHT, BOTTOM, TOP, NEARR, FARR);
    glScalef (2.0, 2.0, 1.0);
    glTranslatef (-1, -1,-(NEARR + (FARR-NEARR)/2)); // glue NEAR to the screen
    // glScalef (1.0, 1.0, 2.0);
    glMatrixMode (GL_MODELVIEW);
    // glLoadIdentity ();
}

static GLvoid glutKeyPressed(unsigned char key, int, int)
{
    static int x = 0, y = 0, w = 320, h = 200;
    if ('p' == key) paused = ! paused;
    else if ('q' == key) exit (0);
    //D else if ('a' == key) ry += 1.f;
    //D else if ('d' == key) ry -= 1.f;
    else if (key == 'f') {
        if (!fs) {
            x = glutGet (GLUT_WINDOW_X), y = glutGet (GLUT_WINDOW_Y);
            w = glutGet (GLUT_WINDOW_WIDTH), h = glutGet (GLUT_WINDOW_HEIGHT);
            glutFullScreen ();
            fs = true;
            mouse_moved = time (0);
        } else
            glutSetCursor (GLUT_CURSOR_INHERIT),
            glutReshapeWindow (w, h), glutPositionWindow (x, y), fs = false;
    }
}

static GLvoid glutPassiveMotion(int a, int b)
{
    static int _pa {}, _pb {};
    if (_pa != a || _pb != b) {
        mouse_moved = time (0);
        glutSetCursor (GLUT_CURSOR_INHERIT);
    }
    _pa = a; _pb = b;
}

// Bresenham's line algorithm. Nothing scary - its just fast.
static inline void Swap(int * a, int * b) { int t = *a; *a = *b; *b = t; }
static inline void Line(int x0, int y0, int x1, int y1, unsigned char * bitmap,
    int w, int h, unsigned int v)
{
    (void)h;
    int step = abs (y1 - y0) > abs (x1 - x0);
    if (step) {
        Swap (&x0, &y0);
        Swap (&x1, &y1);
    }
    if (x0 > x1) {
        Swap (&x0, &x1);
        Swap (&y0, &y1);
    }
    int deltax = x1 - x0;
    int deltay = abs(y1 - y0);
    int error = deltax / 2;
    int ystep;
    int y = y0;
    int x;
    if (y0 < y1) ystep = 1;
    else ystep = -1;
    for (x = x0; x <= x1; x++) {
        if (step && x*w+y >= 0 && x*w+y < w*h) {
            ((unsigned int *)bitmap)[x*w + y] =  v;
            /*for (int i = 1 ; i < 3; i++) {
                if (x * w + y+i < w*h) ((unsigned int *)bitmap)[x*w + y+i] = v;
                if (x * w + y-i >= 0) ((unsigned int *)bitmap)[x*w + y-i] = v;
            }*/
        }
        else if (y*w+x >= 0 && y*w+x < w*h) {
            ((unsigned int *)bitmap)[y*w + x] = v;
            /*for (int i = 1 ; i < 3; i++) {
                if (y * w + x+i < w*h) ((unsigned int *)bitmap)[y*w + x+i] = v;
                if (y * w + x-i >= 0) ((unsigned int *)bitmap)[y*w + x-i] = v;
            }*/
        }
        error = error - deltay;
        if (error < 0) {
            y = y + ystep;
            error = error + deltax;
        }
    }
}// Line()

// Something-like-a-snowflake generator; generates to an RGBA buffer of size "a"
// with alpha = 0xaa.
class SnoflakeGen final
{
    private int _a; // size "a"
    private static int const _I {2}; // max iterations
    struct P final // 2D point/vector - whatever role is required
    {
        int _x, _y;
        P& operator+=(const P & p) { _x += p._x; _y += p._y; return *this; }
        P operator+(const P & p) const { return P {*this} += p; }
        P Translate(int x, int y) const
        {
            P result {*this};
            result._x += x;
            result._y += y;
            return result;
        }
        // Rotate _one above, "a" [degrees], around "p". "p" is expected
        // to be translated to 0,0 prior this call.
        P Rotate(const P & p, float a) const
        {
            // y = a*x + b
            float a_radians = a*MY_PI/180.0;// degrees = radians*180.0/3.14
            float p_tita = p._x ? p._y
//TODO to os.h
#ifdef _WIN32
                ? atan (fabs (static_cast<float>(p._y)/p._x)) : 0.0 : MY_PI/2.0;
#else
                ? atan (abs (static_cast<float>(p._y)/p._x)) : 0.0 : MY_PI/2.0;
#endif
            P result {*this};
            // printf ("p.x: %d, p.y: %d, arctan: %5.2f\n", p._x, p._y, p_tita);
            result._x = static_cast<int>(round (result._x *
                (p._x < 0 ? -1.0 : 1.0) * cos (p_tita-a_radians)));
            result._y = static_cast<int>(round (result._y *
                (p._y < 0 ? -1.0 : 1.0) * sin (p_tita-a_radians)));
            return result;
        }
        P(size_t a) : _x{static_cast<int>(a)}, _y{_x} {}
        // x1 < x2 && y1 < y2
        bool Outside(int x1, int y1, int x2, int y2) const
        {
            return _x < x1 || _x > x2 || _y < y1 || _y > y2;
        }
        P & operator=(const P & p) { _x = p._x; _y = p._y; return *this; }
        P(const P & p) { this->operator= (p); }
        P(int x, int y) : _x{x}, _y{y} {}
        const int & X() { return _x; }
        const int & Y() { return _y; }
    };// class P final
    private P _one, _c; // _one: line segment length; _c: center
    //TODO there is no reason "_one" to be a "P"
    private DMem<unsigned char> _bmp; // RGBA back buffer
    private unsigned int _color {0xaaaaaaffu}; // guess
    private void Line(const P & p1, const P & p2)
    {
        ::Line (p1._x, p1._y, p2._x, p2._y, _bmp, _a, _a, _color);
    }
    private void Gen()
    {
        //D border
        //D ::Line (0, 0, 0, _a-1, _bmp, _a, _a, 0xff808080u);
        //D ::Line (0, _a-1, _a-1, _a-1, _bmp, _a, _a, 0xff808080u);
        //D ::Line (_a-1, _a-1, _a-1, 0, _bmp, _a, _a, 0xff808080u);
        //D ::Line (_a-1, 0, 0, 0, _bmp, _a, _a, 0xff808080u);

        P up_vector {0, _c._y}; // 12 o'clock
        // dance
        for (int i = 0; i < 360; i += 60) {
            P n {_c + _one.Rotate (up_vector, i)};
            Line (_c, n);
            Gen (n, _c, 0);
        }
    }
    // p: newest point; pp: prev point; l: recursion depth
    private void Gen(P & p, P & pp, int l)
    {
        if (p.Outside (0, 0, _a-1, _a-1)) return;
        if (_I == l) return;
        int d = 0;// random () % 2; perhaps someday, but not now
        switch (d) {
            // case 0: return;
            case 0: { // extend by _one
                int tx = p._x, ty = p._y;
                int one_x = p._x == pp._x ? 0
                    : p._x < pp._x ? -_one._x : _one._x,
                    one_y = p._y == pp._y ? 0
                    : p._y < pp._y ? -_one._y : _one._y;
                int nx = abs(p._x - pp._x), ny = abs(p._y - pp._y);
                if (nx > ny) { // y has partial pixels
                   tx += one_x;
                   ty += nx ? static_cast<int>(round (
                       static_cast<float>(ny)/nx * one_y)) : 0;
                }
                else {
                   tx += ny ? static_cast<int>(round (
                       static_cast<float>(nx)/ny * one_x)) : 0;
                   ty += one_y;
                }
                P p2 {tx, ty};
                Line (p, p2);
                Gen (p2, p, l+1);
                //D if (l > 0) return;
            };
            case 1: { // branch to 2 lines at +60 and -60 degrees
                auto p2 = p + _one.Rotate (p.Translate (-_c._x, -_c._y),  60);
                auto p3 = p + _one.Rotate (p.Translate (-_c._x, -_c._y), -60);
                Line (p, p2); Line (p, p3);
                Gen (p2, p, l+1);
                Gen (p3, p, l+1);
                return;
            } break;
            default: return;
        }
    }// Gen(p, pp, l)
    public SnoflakeGen(size_t a)
        : _a{static_cast<int>(a)}, _one {a/6}, _c {a/2}, _bmp {a*a*4}
    {
        Gen ();
    }
    public const unsigned char * Bmp() { return _bmp; }
};// SnowflakeGen

namespace __pointless_verbosity
{
    struct try_finally_gl_res__ final
    {
        GLsizei _n;
        GLuint * _t;
        try_finally_gl_res__(GLsizei n, GLuint * t) : _n{n}, _t{t} {}
        ~try_finally_gl_res__() { if (_t && _n) glDeleteTextures (_n, _t); }
    };
}

int main()
{
    /*float m {}, n {}; //TODO apply the noise to wx, wy, wz above
    for (int i = 0; i < WT_NUM; i++)
        WTable[i] = Noise(i, WT_NUM, m, n);
    for (int i = 0; i < WT_NUM; i++) {
        WTable[i] = (WTable[i] - n)  / (m - n); // this is part of the noise
        printf ("%5.2f ", WTable[i]);
        if (not ((i+1) % WT_SIZE)) printf ("\n");
    }
    {
        DMem<unsigned char> pixel_buf_store {WT_NUM<<2};
        unsigned char * pixel_buf = pixel_buf_store;
        for (int  i = 0; i < WT_NUM; i++) {
            pixel_buf[i*4+0] = pixel_buf[i*4+1] = pixel_buf[i*4+2] =
                static_cast<int>(roundf (WTable[i]*255));
            pixel_buf[i*4+3] = 255;
        }
        BMP (WT_SIZE, WT_SIZE, pixel_buf, "noise.bmp");
    }*/

    const int TEX0_SZ {128};
    SnoflakeGen tex0 {TEX0_SZ};
    // BMP(TEX0_SZ, TEX0_SZ, tex0.Bmp (), "tex0.bmp");

    int GLUT_ARGC = 0;
    glutInit (&GLUT_ARGC, 0);
    glutSetOption (GLUT_ACTION_ON_WINDOW_CLOSE,
                   GLUT_ACTION_GLUTMAINLOOP_RETURNS);
    current_viewport_width = glutGet (GLUT_SCREEN_WIDTH);
    current_viewport_height = glutGet (GLUT_SCREEN_HEIGHT);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_ALPHA);
    glutInitWindowSize (current_viewport_width, current_viewport_height);
    glutInitWindowPosition (0, 0);
    glutCreateWindow ("snow gl 1.0");
    // glEnable (GL_COLOR_MATERIAL);

    glEnable (GL_TEXTURE_2D);

    GLuint t;
    glGenTextures (1, &t);
    __pointless_verbosity::try_finally_gl_res__ ____ {1, &t}; // not Morse code
    glBindTexture (GL_TEXTURE_2D, t);
    float tex_filtering = 0.0f;
    glGetFloatv (GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &tex_filtering);
    // printf ("Anisotropy filter max: %5.2f\n", tex_filtering);
    glTexParameterf (
        GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, tex_filtering);

    //TODO bell-blur the tex_t.Bmp (); say 0.1 0.1 0.2 0.3 0.5 0.8
    glTexParameteri (
        GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexEnvi (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    for (int i = TEX0_SZ, j = 0; i > 0; i >>= 1) {
        SnoflakeGen tex_t {static_cast<size_t>(i)};
        glTexImage2D (GL_TEXTURE_2D, j++, GL_RGBA, i, i, 0, GL_BGRA,
            GL_UNSIGNED_BYTE, tex_t.Bmp ());
    }

    glDisable (GL_LIGHTING);
    glShadeModel (GL_FLAT);
    glClearColor (0.0f, 0.0f, 0.0f, 0.0f);
    glClearDepth (1.0f);
    glDisable (GL_DEPTH_TEST);
    glEnable (GL_POINT_SMOOTH);
    glEnable (GL_ALPHA_TEST);
    glAlphaFunc (GL_GEQUAL, 0.01f);
    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glutKeyboardFunc (glutKeyPressed);
    glutReshapeFunc (glutReshape);
    glutDisplayFunc (glutDisplay);
    glutIdleFunc (glutIdle);
    glutPassiveMotionFunc (glutPassiveMotion);
    mouse_moved = time (nullptr);
    glutMainLoop ();
    return 0;
}