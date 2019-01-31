//---------------------------------------------------------------------------
//algebraic surfaces in parametric form
//---------------------------------------------------------------------------

package com.vs.algebraicsurfaces.algebraicsurfaces;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.util.ArrayList;

import javax.microedition.khronos.opengles.GL10;

public abstract class ParametricSurface {

    int res = 35, nCoords;
    float fromU, toU, fromV, toV, difU, difV;
    float mx, my, mz, mix, miy, miz, MaxVal, MinVal, Dif;
    boolean scaled = false;
    float MAXfloat = Float.MAX_VALUE, Pi = (float)Math.PI, twoPi = Pi * 2f;
    XYZ p = new XYZ();
    ArrayList<Float> lTextures = new ArrayList<Float>(); // textures & vertex list
    ArrayList<XYZ> lVertex = new ArrayList<XYZ>();
    FloatBuffer fptx, fpvx;

    public float sqr(float x) {
        return x * x;
    }

    public float sqr3(float x) {
        return x * x * x;
    }

    public float cube(float x) {
        return x * x * x;
    }

    public float sqr4(float x) {
        return x * x * x * x;
    }

    public float sqr5(float x) {
        return x * x * x * x * x;
    }

    public float fabs(float x) {
        return (x < 0) ? -x : x;
    }

    public float sin(float x) {
        return (float) Math.sin(x);
    }

    public float cos(float x) {
        return (float) Math.cos(x);
    }

    public float exp(float x) {
        return (float) Math.exp(x);
    }

    public float log(float x) {
        return (float) Math.log(x);
    }

    public float pow(float x, float y) {
        return (float) Math.pow(x, y);
    }

    public float sinh(float x) {
        return (float) Math.sinh(x);
    }

    public float cosh(float x) {
        return (float) Math.cosh(x);
    }

    ////////////////////////////////////////////////////////////////////////////////
    //Base Parametric Surface class
    //	 required implementing virtual func: XYZ ParametricFunc(float u, float v)
    //	 p.x=x(u,v)
    //	 p.y=y(u,v)
    //	 p.z=z(u,v)
    ////////////////////////////////////////////////////////////////////////////////
    // default func returns a 0 point
    abstract XYZ ParametricFunc(float u, float v); // assign 'p' with x,y,z

    public void setResol(int resol) {
        res = resol;
    }

    public float rad2Deg(float rad) {
        return rad * 180f / Pi;
    }

    private float scaleU(float val) {
        return val * difU + fromU;
    }

    private float scaleV(float val) {
        return val * difV + fromV;
    }

    private float scale01(float val) {
        return val / Dif;
    } // keep center

    private float scale0to1(float val) {
        return (val - MinVal) / Dif;
    } // scale 0..1

    private void scale() {
        p.x = scale01(p.x);
        p.y = scale01(p.y);
        p.z = scale01(p.z);
    }

    XYZ Eval(float u, float v) {
        p = ParametricFunc(scaleU(u), scaleV(v));
        return p;
    }

    void minMaxP() { // update min/max in p
        if (p.x > mx) mx = p.x;
        if (p.y > my) my = p.y;
        if (p.z > mz) mz = p.z;
        if (p.x < mix) mix = p.x;
        if (p.y < miy) miy = p.y;
        if (p.z < miz) miz = p.z;
    }

    void calcDif() {
        if (mx < my) {
            if (mx < mz) MaxVal = mx;
            else MaxVal = mz;
        } else {
            if (my < mz) MaxVal = my;
            else MaxVal = mz;
        }
        if (mix < miy) {
            if (mix < miz) MinVal = mix;
            else MinVal = miz;
        } else {
            if (miy < miz) MinVal = miy;
            else MinVal = miz;
        }

        Dif = fabs(MaxVal - MinVal);
    }

    void glTexCoord2f(float tx, float ty) {
        lTextures.add(tx);
        lTextures.add(ty);
    } // add textures & vertex coords

    void glVertex3f(XYZ p) {
        lVertex.add(new XYZ(p));
    }

    void addTextVertex(float tx, float ty) {  // add textures coords and vertex
        glTexCoord2f(tx, ty);
        glVertex3f(Eval(tx, ty));
        minMaxP();
    }

    void scaleCoords() {
        nCoords = lVertex.size();
        calcDif();
        for (int i = 0; i < nCoords; i++) {
            p = lVertex.get(i);
            scale();
            lVertex.set(i, p);
        }
    }

    void initMinMax() {
        mx = my = mz = -MAXfloat;
        mix = miy = miz = MAXfloat;
    }

    void initLists() {
        lTextures.clear();
        lVertex.clear();
        fptx = fpvx = null;
    }

    void calcCoords(int resol, float _fromU, float _toU, float _fromV, float _toV) { // calc & load
        if (res == resol) return; // calculated? -> yes:return

        fromU = _fromU;
        toU = _toU;
        fromV = _fromV;
        toV = _toV; // define limits
        difU = fabs(fromU - toU);
        difV = fabs(fromV - toV);

        setResol(resol);
        initLists();
        initMinMax();

        float dr = 1f / res, dt = 1f / res;
        for (int i = 0; i < res; i++) { // generated QUADS
            float idr = i * dr;
            for (int j = 0; j < res; j++) {
                float jdt = j * dt, jdr = jdt;
                addTextVertex(idr, jdr);
                addTextVertex(idr + dr, jdt);
                addTextVertex(idr + dr, jdt + dt);
                addTextVertex(idr, jdt + dt);
            }
        }
        scaleCoords();
        loadCoords(); // fp <- list
    }

    void loadCoords() {
        float[] tx = new float[lTextures.size()];
        for (int i = 0; i < lTextures.size(); i++) tx[i] = lTextures.get(i);
        nCoords = lVertex.size();
        float[] vx = new float[nCoords * 3];
        for (int i = 0; i < nCoords; i++) {
            p = lVertex.get(i);
            vx[i * 3 + 0] = p.x;
            vx[i * 3 + 1] = p.y;
            vx[i * 3 + 2] = p.z;
        }
        fptx = putCoords(tx);
        fpvx = putCoords(vx);
    }

    FloatBuffer putCoords(float[] coords) {
        FloatBuffer pntBuff = (ByteBuffer.allocateDirect(coords.length * 4).order(ByteOrder.nativeOrder())).asFloatBuffer();
        pntBuff.put(coords).position(0);
        return pntBuff;
    }

    void draw(GL10 gl, boolean textures) {
        if (fpvx == null || fptx == null) return; // calculated coords??

        gl.glVertexPointer(3, GL10.GL_FLOAT, 0, fpvx);
        if (textures) gl.glTexCoordPointer(2, GL10.GL_FLOAT, 0, fptx);
        for (int i = 0; i < nCoords; i += 4) { // draw QUADS of coords
            gl.glDrawArrays(GL10.GL_LINE_LOOP, i, 4);
            if (textures) gl.glDrawArrays(GL10.GL_TRIANGLE_FAN, i, 4);
        }
    }

    class XYZ {
        float x, y, z;

        XYZ() {
            x = y = z = 0;
        }

        XYZ(XYZ p) {
            this.x = p.x;
            this.y = p.y;
            this.z = p.z;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
//surface implementation
////////////////////////////////////////////////////////////////////////////////

class Ttanaka extends ParametricSurface {
    float[][] paramSets = {{0, 4, 3, 4, 5, 7, 4}, {0, 4, 3, 0, 5, 7, 4}, {0, 3, 4, 8, 5, 5, 2}, {14, 3, 1, 8, 5, 5, 2}};
    float a = 0, // center hole size of a torus
            b1 = 4,// number of cross
            b2 = 3,// number of cross
            c = 4, // distance from the center of rotation
            d = 5, // number of torus
            w = 7, // gap width
            h = 4; // height

    void setParam(int np) {
        a = paramSets[np][0];
        b1 = paramSets[np][1];
        b2 = paramSets[np][2];
        c = paramSets[np][3];
        d = paramSets[np][4];
        w = paramSets[np][5];
        h = paramSets[np][6];
    }

    float getNTorus() {
        return d;
    } // number of torus

    float f(float v) {
        return sin(2 * sin(sin(sin(v))));
    }

    XYZ ParametricFunc(float s, float t) {
        p.x = (a - cos(t) + w * sin(b1 * s)) * cos(b2 * s);
        p.y = (a - cos(t) + w * sin(b1 * s)) * f(b2 * s);
        p.z = h * (w * sin(b1 * s) + f(t)) + c;
        return p;
    }
}

class TKuen extends ParametricSurface {
    XYZ ParametricFunc(float u, float v) {
        p.x = 2 * cosh(v) * (cos(u) + u * sin(u)) / (cosh(v) * cosh(v) + u * u);
        p.y = 2 * cosh(v) * (-u * cos(u) + sin(u)) / (cosh(v) * cosh(v) + u * u);
        p.z = v - (2 * sinh(v) * cosh(v)) / (cosh(v) * cosh(v) + u * u);
        return p;
    }
}

class TButterFly extends ParametricSurface {
    XYZ ParametricFunc(float u, float v) {
        float t1 = (exp(cos(u)) - 2 * cos(4 * u) + sqr5(sin(u / 12))) * sin(v);
        p.x = sin(u) * t1;
        p.y = cos(u) * t1;
        p.z = sin(v);
        return p;
    }
}

class TRoseSurface extends ParametricSurface {
    XYZ ParametricFunc(float u, float v) {
        float a = 1, n = 7;
        p.x = a * sin(n * u) * cos(u) * sin(v);
        p.y = a * sin(n * u) * sin(u) * sin(v);
        p.z = cos(v) / (n * 3);
        return p;
    }
}

class TUpDownShellSurface extends ParametricSurface {
    XYZ ParametricFunc(float u, float v) {
        p.x = u * sin(u) * cos(v);
        p.y = u * cos(u) * cos(v);
        p.z = u * sin(v); // -10,10, -10,10
        return p;
    }
}


class TCayleySurface extends ParametricSurface {
    XYZ ParametricFunc(float u, float v) {
        // plot3d([u*sin(v)-u*cos(v), u^2*sin(v)*cos(v), u^3*sin(v)^2*cos(v)], u=0..0.5, v=0..2*Pi, numpoints=1000);
        p.x = u * sin(v) - u * cos(v);
        p.y = sqr(u) * sin(v) * cos(v);
        p.z = cube(u) * sqr(sin(v)) * cos(v);
        return p;
    }
}


class TPluckerConoidSurface extends ParametricSurface {
    XYZ ParametricFunc(float u, float v) {
        p.x = u * v;
        p.y = u * (float) Math.sqrt(1 - sqr(v));
        p.z = 1 - sqr(v);
        return p;
    }
}


class TAmmoniteSurface extends ParametricSurface {
    float W(float u) {
        return pow(u / (2 * Pi), 2.2f);
    }

    XYZ ParametricFunc(float u, float v) {
        float N = 5.6f;  // number of turns
        float F = 120f;  // wave frequency
        float A = 0.2f;  // wave amplitude

        p.x = W(u) * cos(N * u) * (2 + sin(v + cos(F * u) * A));
        p.y = W(u) * sin(N * u) * (2 + sin(v + cos(F * u) * A));
        p.z = W(u) * cos(v);
        return p;
    }
}

class TAppleSurface extends ParametricSurface {
    XYZ ParametricFunc(float u, float v) {
        float R1 = 4, R2 = 3.8f;

        p.x = cos(u) * (R1 + R2 * cos(v)) + pow((v / Pi), 100);
        p.y = sin(u) * (R1 + R2 * cos(v)) + 0.25f * cos(5 * u);
        p.z = -2.3f * log(1 - v * 0.3157f) + 6 * sin(v) + 2 * cos(v);
        return p;
    }
}


class TAstroidalEllipseSurface extends ParametricSurface {
    XYZ ParametricFunc(float u, float v) {
        float A = 1, B = 1, C = 1;

        p.x = pow(A * cos(u) * cos(v), 3);
        p.y = pow(B * sin(u) * cos(v), 3);
        p.z = pow(C * sin(v), 3);
        return p;
    }
    // <0,0>,<2*pi,2*pi>
}

class TBohemianDomeSurface extends ParametricSurface {
    XYZ ParametricFunc(float u, float v) {
        float A = 0.5f, B = 1.5f, C = 1;
        p.x = A * cos(u);
        p.y = B * cos(v) + A * sin(u);
        p.z = C * sin(v);
        return p;
    }
}

class TConicalSpiralSurface extends ParametricSurface {
    XYZ ParametricFunc(float u, float v) {
        p.x = u * v * sin(15 * v);
        p.y = v;
        p.z = u * v * cos(15 * v);
        return p;
        // <0,-1>,<1,1>
    }
}

class TEnneperSurface extends ParametricSurface {
    XYZ ParametricFunc(float u, float v) {
        p.x = u - u * u * u / 3 + u * v * v;
        p.y = v - v * v * v / 3 + v * u * u;
        p.z = u * u - v * v;
        return p;
    }
}


class TScherk extends ParametricSurface {
    XYZ ParametricFunc(float u, float v) {
        float aa = 0.1f;
        v += 0.1;
        p.x = u;
        p.y = v;
        p.z = (log(fabs(cos(aa * v) / cos(aa * u)))) / aa;
        return p;
    }
}


class TDini extends ParametricSurface {
    XYZ ParametricFunc(float u, float v) {
        float psi = 0.3f; // aa;
        if (psi < 0.001f) psi = 0.001f;
        if (psi > 0.999f) psi = 0.999f;
        psi = psi * Pi;
        float sinpsi = sin(psi);
        float cospsi = cos(psi);
        float g = (u - cospsi * v) / sinpsi;
        float s = exp(g);
        float r = (2 * sinpsi) / (s + 1 / s);
        float t = r * (s - 1 / s) * 0.5f;

        p.x = (u - t);
        p.y = (r * cos(v));
        p.z = (r * sin(v));

        return p;
    }
}

class TKleinBottle extends ParametricSurface {
    float t;

    TKleinBottle() {
        t = 4.5f;
    }

    XYZ ParametricFunc(float u, float v) {
        float tmp = (4 + 2 * cos(u) * cos(t * v) - sin(2 * u) * sin(t * v));
        p.x = sin(v) * tmp;
        p.y = cos(v) * tmp;
        p.z = 2 * cos(u) * sin(t * v) + sin(2 * u) * cos(t * v);
        return p;
    }
}

class TKleinBottle0 extends ParametricSurface {
    float t;

    TKleinBottle0() {
        t = 4.5f;
    }

    XYZ ParametricFunc(float u, float v) {
        p.x = (0 <= u && u < Pi) ? 6 * cos(u) * (1 + sin(u)) + 4 * (1 - 0.5f * cos(u)) * cos(u) * cos(v) : 6 * cos(u) * (1 + sin(u)) + 4 * (1 - 0.5f * cos(u)) * cos(v + Pi);
        p.y = (0 <= u && u < Pi) ? 16 * sin(u) + 4 * (1 - 0.5f * cos(u)) * sin(u) * cos(v) : 16 * sin(u);
        p.z = 4 * (1 - 0.5f * cos(u)) * sin(v);
        return p;
    }
}

class TBour extends ParametricSurface {
    XYZ ParametricFunc(float u, float v) {
        p.x = u * cos(v) - 0.5f * u * u * cos(2 * v);
        p.y = -u * sin(v) - 0.5f * u * u * sin(2 * v);
        p.z = 4 / 3 * pow(u, 1.5f) * cos(1.5f * v);
        return p;
    }
}

class TBreatherSurface extends ParametricSurface {
    float aa, w1, w;

    TBreatherSurface() {
        aa = 0.45f;  // Values from 0.4 to 0.6 produce sensible results
        w1 = 1 - aa * aa;
        w = (float) Math.sqrt(w1);
    }

    float d(float u, float v) {
        return aa * (pow((w * cosh(aa * u)), 2) + pow((aa * sin(w * v)), 2));
    }

    XYZ ParametricFunc(float u, float v) {
        p.x = -u + (2 * w1 * cosh(aa * u) * sinh(aa * u) / d(u, v));
        p.y = 2 * w * cosh(aa * u) * (-(w * cos(v) * cos(w * v)) - (sin(v) * sin(w * v))) / d(u, v);
        p.z = 2 * w * cosh(aa * u) * (-(w * sin(v) * cos(w * v)) + (cos(v) * sin(w * v))) / d(u, v);
        return p;
    }

}


/*
Seashell surfaces
These are parametric isosurfaces using variations of functions like:

#declare N=5.6;  // number of turns
#declare H=3.5;  // height
#declare P=2;    // power
#declare L=4;    // Controls spike length
#declare K=9;    // Controls spike sharpness
 */

class TSeashell extends ParametricSurface {
    float N, H, P, L, K;

    TSeashell() {
        N = 5.6f;  // number of turns
        H = 3.5f;  // height
        P = 2;    // power
        L = 4;    // Controls spike length
        K = 9;
    }

    float W(float u) {
        return pow(u / (2 * Pi), P);
    }

    XYZ ParametricFunc(float u, float v) {
        p.x = W(u) * cos(N * u) * (1 + cos(v));
        p.y = W(u) * sin(N * u) * (1 + cos(v));
        p.z = W(u) * (sin(v) + pow(sin(v / 2), K) * L) + H * pow(u / (2 * Pi), P + 1);
        return p;
    }
}

//TudoRose surface
class TTudorRoseSurface extends ParametricSurface {

    float R(float u, float v) {
        return cos(v) * cos(v) * Math.max(fabs(sin(4 * u)), 0.9f - 0.2f * fabs(cos(8 * u)));
    }

    XYZ ParametricFunc(float u, float v) {
        p.x = R(u, v) * cos(u) * cos(v);
        p.y = R(u, v) * sin(u) * cos(v);
        p.z = R(u, v) * sin(v) * 0.5f;
        return p;
    }
}

//Cap's surface
class TCapSurface extends ParametricSurface {
    XYZ ParametricFunc(float u, float v) {
        p.x = 0.5f * cos(u) * sin(2 * v);
        p.y = 0.5f * sin(u) * sin(2 * v);
        p.z = 0.5f * (sqr(cos(v)) - sqr(cos(u)) * sqr(sin(v)));
        return p;
    }
}

//Boy's surface
class TBoySurface extends ParametricSurface {
    XYZ ParametricFunc(float u, float v) {
        float dv = (2 - (float) Math.sqrt(2) * sin(3 * u) * sin(2 * v)),
                d1 = (cos(u) * sin(2 * v)),
                d2 = (float) Math.sqrt(2) * sqr(cos(v));

        p.x = (d2 * cos(2 * u) + d1) / dv;
        p.y = (d2 * sin(2 * u) + d1) / dv;
        p.z = (3 * sqr(cos(v))) / (2 - (float) Math.sqrt(2) * sin(3 * u) * sin(2 * v));

		/* x,y,z deformations
		float a=2.1,b=1.3,c=1;
		p.x*=a; p.y*=b; p.z*=c;
		 */
        return p;
    }
}

//
//roman algebraic surface
//
class TRomanSurface extends ParametricSurface {
    XYZ ParametricFunc(float r, float t) // 0 <= r <= 1   |  0 <= t <= 2*PI
    {
        float r2 = r * r, rq = (float) Math.sqrt((1 - r2)), st = sin(t), ct = cos(t);
        p.x = r2 * st * ct;
        p.y = r * st * rq;
        p.z = r * ct * rq;
        return p;
    }
}


class SurfaceTest { // TODO: add here new funcs
    //----------------------------------------------------------------------------------------------------------------------------
    public String srfNames[] = {"cap", "boy", "roman", "sea shell", "tudor rose", "breather", "klein bottle",
            "klein bottle 0", "bour", "dini", "enneper", "scherk", "conical spiral", "bohemian dome", "astrodial ellipse",
            "apple", "ammonite", "plucker comoid", "cayley", "up down shell", "butterfly", "rose", "kuen"
            , "tanaka-0", "tanaka-1", "tanaka-2", "tanaka-3"};
    float Pi = 3.141592f, twoPi = Pi * 2;
    //----------------------------------------------------------------------------------------------------------------------------
    TRoseSurface RoseSrf = new TRoseSurface();
    TRomanSurface Rsrf = new TRomanSurface();
    TCapSurface Csrf = new TCapSurface();
    TBoySurface Bsrf = new TBoySurface();
    TSeashell Ssrf = new TSeashell();
    TTudorRoseSurface Tsrf = new TTudorRoseSurface();
    TBreatherSurface Brsrf = new TBreatherSurface();
    TKleinBottle KBsrf = new TKleinBottle();
    TKleinBottle0 KBsrf0 = new TKleinBottle0();
    TBour BourS = new TBour();
    TDini DiniS = new TDini();
    TEnneperSurface EnneperS = new TEnneperSurface();
    TScherk Scherk = new TScherk();
    TConicalSpiralSurface ConSpi = new TConicalSpiralSurface();
    TBohemianDomeSurface BDsurf = new TBohemianDomeSurface();
    TAstroidalEllipseSurface AESurf = new TAstroidalEllipseSurface();
    TAppleSurface SApple = new TAppleSurface();
    TAmmoniteSurface AmmSrf = new TAmmoniteSurface();
    TPluckerConoidSurface PCSrf = new TPluckerConoidSurface();
    TCayleySurface CYSrf = new TCayleySurface();
    TUpDownShellSurface UDSSrf = new TUpDownShellSurface();
    TButterFly bfy = new TButterFly();
    TKuen KuenSrf = new TKuen();
    Ttanaka tanaka = new Ttanaka();
    ParametricSurface[] ps = {Csrf, Rsrf, RoseSrf};

    public void calcCoords(int ns, int resol) {
        switch (ns) {
            case 0:
                Csrf.calcCoords(resol, 0, Pi, 0, Pi);
                break;
            case 1:
                Bsrf.calcCoords(resol, 0, Pi, 0, Pi);
                break;
            case 2:
                Rsrf.calcCoords(resol, 0, 1, 0, twoPi);
                break;
            case 3:
                Ssrf.calcCoords(resol, 0, twoPi, 0, twoPi);
                break;
            case 4:
                Tsrf.calcCoords(resol, 0, Pi, 0, Pi);
                break;
            case 5:
                Brsrf.calcCoords(resol, -20, 20, 20, 80);
                break;
            case 6:
                KBsrf.calcCoords(resol, 0, twoPi, 0, twoPi);
                break;
            case 7:
                KBsrf0.calcCoords(resol, 0, twoPi, 0, twoPi);
                break;
            case 8:
                BourS.calcCoords(resol, 0, twoPi, 0, twoPi);
                break;
            case 9:
                DiniS.calcCoords(resol, 0, twoPi, 0, twoPi);
                break;
            case 10:
                EnneperS.calcCoords(resol, -1, 1, -1, 1);
                break;
            case 11:
                Scherk.calcCoords(resol, 1, 30, 1, 30);
                break;
            case 12:
                ConSpi.calcCoords(resol, 0, 1, -1, 1);
                break;
            case 13:
                BDsurf.calcCoords(resol, 0, twoPi, 0, twoPi);
                break;
            case 14:
                AESurf.calcCoords(resol, 0, twoPi, 0, twoPi);
                break;
            case 15:
                SApple.calcCoords(resol, 0, twoPi, -Pi, Pi);
                break;
            case 16:
                AmmSrf.calcCoords(resol, 0, twoPi, 0, twoPi);
                break;
            case 17:
                PCSrf.calcCoords(resol, -2, 2, -1, 1);
                break;
            case 18:
                CYSrf.calcCoords(resol, 0, 3, 0, twoPi);
                break;
            case 19:
                UDSSrf.calcCoords(resol, -10, 10, -10, 10);
                break;
            case 20:
                bfy.calcCoords(resol, 0, twoPi, 0, twoPi);
                break;
            case 21:
                RoseSrf.calcCoords(resol, 0, twoPi, 0, twoPi);
                break;
            case 22:
                KuenSrf.calcCoords(resol, -4, 4, -3.75f, +3.75f);
                break;
            case 23:
            case 24:
            case 25:
            case 26:
                tanaka = new Ttanaka();
                tanaka.setParam(ns - 23);
                tanaka.calcCoords(resol, 0, twoPi, 0, twoPi);
                break;
        }
    }

    public void draw(GL10 gl, int ns, boolean textures) {
        switch (ns) {
            case 0:
                Csrf.draw(gl, textures);
                break;
            case 1:
                Bsrf.draw(gl, textures);
                break;
            case 2:
                Rsrf.draw(gl, textures);
                break;
            case 3:
                Ssrf.draw(gl, textures);
                break;
            case 4:
                Tsrf.draw(gl, textures);
                break;
            case 5:
                Brsrf.draw(gl, textures);
                break;
            case 6:
                KBsrf.draw(gl, textures);
                break;
            case 7:
                KBsrf0.draw(gl, textures);
                break;
            case 8:
                BourS.draw(gl, textures);
                break;
            case 9:
                DiniS.draw(gl, textures);
                break;
            case 10:
                EnneperS.draw(gl, textures);
                break;
            case 11:
                Scherk.draw(gl, textures);
                break;
            case 12:
                ConSpi.draw(gl, textures);
                break;
            case 13:
                BDsurf.draw(gl, textures);
                break;
            case 14:
                AESurf.draw(gl, textures);
                break;
            case 15:
                SApple.draw(gl, textures);
                break;
            case 16:
                AmmSrf.draw(gl, textures);
                break;
            case 17:
                PCSrf.draw(gl, textures);
                break;
            case 18:
                CYSrf.draw(gl, textures);
                break;
            case 19:
                UDSSrf.draw(gl, textures);
                break;
            case 20:
                bfy.draw(gl, textures);
                break;
            case 21:
                RoseSrf.draw(gl, textures);
                break;
            case 22:
                KuenSrf.draw(gl, textures);
                break;
            case 23:
            case 24:
            case 25:
            case 26:
                for (int i = 0; i < tanaka.getNTorus(); i++) { // rm = Table[{x, y, z}.RotationMatrix[2 i Pi/d, {1, 1, 1}], {i, d}];
                    tanaka.draw(gl, textures);                // ParametricPlot3D[rm, {t, 0, 2 Pi}, {s, 0, 2 Pi}]
                    gl.glRotatef(tanaka.rad2Deg(twoPi / tanaka.getNTorus()), 1, 1, 1);
                }
                break;
        }
    }
}
