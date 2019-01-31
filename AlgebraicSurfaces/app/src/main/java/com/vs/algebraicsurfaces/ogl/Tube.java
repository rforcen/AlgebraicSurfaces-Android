package com.vs.algebraicsurfaces.ogl;

import javax.microedition.khronos.opengles.GL10;

public class Tube {

    Vertex root, end;
    Quad[] polys;
    float r;
    int slices, nVertex;
    fbVertex fbCoords, fbTextures, fbColors, fbNormals;
    float[] defaultColor;

    public Tube() {
    }

    public Tube(float[] _root, float[] _end, float r) {
        dotube(new Vertex(_root), new Vertex(_end), r, 8, null, 1.1f);
    }

    public Tube(float[] _root, float[] _end, float r, int slices) {
        dotube(new Vertex(_root), new Vertex(_end), r, slices, null, 1.1f);
    }

    public Tube(float[] _root, float[] _end, float r, int slices, float[] color, float extraLen) {
        dotube(new Vertex(_root), new Vertex(_end), r, slices, color, extraLen);
    }

    public void dotube(Vertex _root, Vertex _end, float _r, int slices, float[] color, float extraLen) {
        polys = new Quad[slices + 1]; // complete with '0' slice
        root = _root;
        end = _end;
        r = _r;
        this.slices = slices;
        nVertex = (slices + 1) * 4;
        defaultColor = color;
        float PI = 3.141592f, HALF_PI = PI / 2f, PIs = PI / (slices / 2); // len*extraLen for bigger overlapping (1.1)
        float len = extraLen * (float) Math.sqrt((end.x - root.x) * (end.x - root.x) + (end.y - root.y) * (end.y - root.y) + (end.z - root.z) * (end.z - root.z));
        len /= 2;
        Vertex vec = new Vertex((end.x - root.x), (end.y - root.y), (end.z - root.z)),
                mid = new Vertex((root.x + end.x) / 2, (root.y + end.y) / 2, ((root.z + end.z) / 2));
        float heading = (float) Math.atan2(vec.z, vec.x);
        heading = 0 - heading;
        float pitch = (float) Math.atan2(vec.y, (float) Math.sqrt((vec.x * vec.x) + (vec.z * vec.z)));

        for (int j = 0; j <= slices; j++) {
            Vertex[] a = new Vertex[4]; // rod is horizontal along X to start
            a[0] = new Vertex(len, r * cos(j * PIs), r * sin(j * PIs));
            a[1] = new Vertex(-len, r * cos(j * PIs), r * sin(j * PIs));
            a[2] = new Vertex(-len, r * cos((j + 1) * PIs), r * sin((j + 1) * PIs));
            a[3] = new Vertex(len, r * cos((j + 1) * PIs), r * sin((j + 1) * PIs));

            for (int i = 0; i < a.length; i++) {//pitch, then heading
                float _x = a[i].x, _y = a[i].y;
                a[i].x = _x * cos(pitch) + _y * cos(pitch + HALF_PI);
                a[i].y = _x * sin(pitch) + _y * sin(pitch + HALF_PI);
            }
            for (int i = 0; i < a.length; i++) {
                float _z = a[i].z, _x = a[i].x;
                a[i].z = _z * cos(heading) + _x * cos(heading + HALF_PI);
                a[i].x = _z * sin(heading) + _x * sin(heading + HALF_PI);
            }
            for (int i = 0; i < a.length; i++) {//move tube so that it goes from one end to the other.
                a[i].x += mid.x;
                a[i].y += mid.y;
                a[i].z += mid.z;
            }
            polys[j] = new Quad(a); // a[0],a[1],a[2],a[3]);
        }
        genVertexBuffers();
    }

    private void genVertexBuffers() {
        fbCoords = new fbVertex(genCoords());
        fbTextures = new fbVertex(genTextures());
        fbColors = new fbVertex(genColors());
        fbNormals = new fbVertex(genNormals());
    }

    private void loadfb(GL10 gl) {
        gl.glVertexPointer(3, GL10.GL_FLOAT, 0, fbCoords.getBuffer());
        gl.glTexCoordPointer(2, GL10.GL_FLOAT, 0, fbTextures.getBuffer());
        gl.glColorPointer(4, GL10.GL_FLOAT, 0, fbColors.getBuffer());
        gl.glNormalPointer(GL10.GL_FLOAT, 0, fbNormals.getBuffer());
    }

    public void draw(GL10 gl) {
        loadfb(gl);
        gl.glDrawArrays(GL10.GL_LINE_LOOP, 0, nVertex);
    }

    public void drawSolid(GL10 gl) {
        loadfb(gl);
        for (int i = 0; i < nVertex; i += 4) gl.glDrawArrays(GL10.GL_TRIANGLE_FAN, i, 4);
    }

    float[][] genCoords() {
        float[][] xyz = new float[(slices + 1) * 4][];
        for (int i = 0; i <= slices; i++) {
            for (int j = 0; j < 4; j++) xyz[i * 4 + j] = polys[i].v[j].getfv(); // 1 quad = 4 vertex
        }
        return xyz;
    }

    float[][] genTextures() {
        float[][] textures = new float[slices + 1][];
        for (int i = 0; i <= slices; i++) {
            float dy = (float) i / slices;
            textures[i] = new float[]{0, 0, 1, 0, 1, dy, 0, dy};
        }
        return textures;
    }

    private float[][] genColors() {
        float[][] xyz = new float[(slices + 1) * 4][];
        for (int i = 0; i <= slices; i++) {
            for (int j = 0; j < 4; j++)
                xyz[i * 4 + j] = (defaultColor == null) ? polys[i].v[j].getColor() : defaultColor; // 1 quad = 4 vertex
        }
        return xyz;
    }

    private float[][] genNormals() {
        float[][] normals = new float[(slices + 1) * 4][];
        for (int i = 0; i <= slices; i++)
            for (int j = 0; j < 4; j++) normals[i * 4 + j] = polys[i].getfNormal();
        return normals;
    }

    int[][] genFaces() { // faces are a sequence 0..polys.length
        int[][] face = new int[slices + 1][];
        for (int i = 0; i <= slices; i++)
            face[i] = new int[]{i * 4 + 0, i * 4 + 1, i * 4 + 2, i * 4 + 3};
        return face;
    }

    float[][] genTrack(float[][] _xyz, float r, int slices) {
        int lt = _xyz.length - 1;
        float[][] xyz = new float[lt * slices * 4][];
        for (int i = 0; i < lt; i++) {
            dotube(new Vertex(_xyz[i]), new Vertex(_xyz[i + 1]), r, slices, defaultColor, 1);
            float[][] gc = genCoords();
            for (int j = 0; j < gc.length; j++)
                xyz[i * slices * 4 + j] = gc[j];
        }
        return xyz;
    }

    private float cos(float x) {
        return (float) Math.cos(x);
    }

    private float sin(float x) {
        return (float) Math.sin(x);
    }
}

class Vertex {    //simple holder for 3 values, saves keeping track of 3 times as many values.
    float x, y, z;

    Vertex() {
    }

    Vertex(float _x, float _y, float _z) {
        x = _x;
        y = _y;
        z = _z;
    }

    Vertex(float[] v) {
        x = v[0];
        y = v[1];
        z = v[2];
    }

    float[] getfv() {
        return new float[]{x, y, z};
    }

    public float[] getColor() {
        return new float[]{1, 1, 1, 1};
    }

    public void set(float _x, float _y, float _z) {
        x = _x;
        y = _y;
        z = _z;
    }

    public void mult(float m) {
        x *= m;
        y *= m;
        z *= m;
    }
}

class Quad {//4 sided poly, again, saves us having to keep track of 12 floats (4 corners, 3 values per vertex...)
    Vertex[] v;
    final int nVertex = 4;
    Vertex pa = new Vertex(), pb = new Vertex(), norm = new Vertex();

    Quad(Vertex a, Vertex b, Vertex c, Vertex d) {
        v = new Vertex[nVertex];
        v[0] = a;
        v[1] = b;
        v[2] = c;
        v[3] = d;
        calcNormal();
    }

    Quad(Vertex[] v) {
        this.v = v;
        calcNormal();
    }

    Vertex getNormal() {
        return norm;
    }

    float[] getfNormal() {
        return new float[]{norm.x, norm.y, norm.z};
    }

    private void calcNormal() { // to norm
        Vertex p = v[0], p1 = v[1], p2 = v[2];
        pa.set(p1.x - p.x, p1.y - p.y, p1.z - p.z);
        pb.set(p2.x - p.x, p2.y - p.y, p2.z - p.z);
        norm.set(pa.y * pb.z - pa.z * pb.y, pa.z * pb.x - pa.x * pb.z, pa.x * pb.y - pa.y * pb.x);
        norm = Normalise(norm);
    }

    private Vertex Normalise(Vertex p) { // Normalise a vector
        double length = p.x * p.x + p.y * p.y + p.z * p.z;
        if (length > 0) {
            length = Math.sqrt(length);
            p.x /= length;
            p.y /= length;
            p.z /= length;
        } else {
            p.x = p.y = p.z = 0;
        }
        return p;
    }
}

class TrackTube {
    private Tube[] tube;

    public void genTrack(float[][] xyz, float tubeRadius, int slices, float[] defaultColor) {
        int lt = xyz.length;
        tube = new Tube[lt];
        for (int i = 0; i < lt - 1; i++)
            tube[i] = new Tube(xyz[i], xyz[i + 1], tubeRadius, slices, defaultColor, 1.1f);
        tube[lt - 1] = new Tube(xyz[0], xyz[lt - 1], tubeRadius, slices, defaultColor, 1.1f);
    }

    public void drawSolid(GL10 gl) {
        for (Tube t : tube) t.drawSolid(gl);
    }

    public void drawWire(GL10 gl) {
        for (Tube t : tube) t.draw(gl);
    }
}