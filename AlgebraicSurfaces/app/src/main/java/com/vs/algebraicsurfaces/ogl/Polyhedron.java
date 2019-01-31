package com.vs.algebraicsurfaces.ogl;

import javax.microedition.khronos.opengles.GL10;
import android.graphics.Color;

public class Polyhedron {       // constructor of initially null Polyhedron
	final float pi=(float)Math.PI, pi2=2*pi, phi=1.618f, sqPhi=(float)Math.sqrt(phi);
	int[][]face;
	float[][]xyz, texture, color, normal;    // array of vertex coords.  xyz.length = # of vertices
	fbVertex fbCoords, fbText, fbColor, fbNormal;
	ColorIndex colIx=new ColorIndex(512, Color.BLUE, Color.RED);

	public enum Type  {Tetrahedron, Cube, Octahedron, Icosahedron, Dodecahedron, Prism, Antiprism, Pyramid, Sphere, Tube, Torus, MultiPanel, SlidePanel, Ring};
	Type type;
	String name = "";
	public Polyhedron(String name) { this.name=name; }
	public Polyhedron() 			{ }
	public Polyhedron(Type type, int n)		{ init(type,n,0);	}
	public Polyhedron(Type type, int w, int h)	{ init(type,w,h);	}
	public Polyhedron(Type type)				{ init(type,4,0);	}
	public int 	getNCoords() 				{ return xyz.length;}
	public float[]	getCoord(int i) 			{ return xyz[i]; }

	void init(Type type, int n, int h)	{
		this.type=type;
		switch (type) {
		case Tetrahedron: 	tetrahedron(); 	break;		case Cube:			cube(); 		break;
		case Octahedron:	octahedron();  	break;		case Icosahedron:	icosahedron(); 	break;
		case Dodecahedron:	dodecahedron(); break;
		case Prism:			prism(n); 		break;		case Antiprism:		antiprism(n); 	break;
		case Pyramid:		pyramid(n); 	break;		case Sphere:		sphere(1,30);	break;
		case Tube:			tube(n,1,.5f);	break;
		case MultiPanel:   	multiPanel(n,h); break;		case SlidePanel:	slidePanel(n,h); break;
		}
	}
	void dimxyz(int n) 		{xyz	=new float[n][];}
	void dimnormal(int n) 	{normal	=new float[n][];}
	void dimface(int n)		{face	=new int[n][];}
	void dimtexture(int n)	{texture=new float[n][];}
	void dimcolor(int n) 	{color	=new float[n][];}
	void addVxyz(int i, float[]v) {xyz[i]=v;}

	void addFaces(int[][]faces) 	{ face=faces; }
	void addxyz(float[][]xyzs) 	{ xyz=xyzs;   }
	private float sqr(float x) {return x*x;}
	// poly operators
	int[]sequence(int start, int stop) {    // make list of integers, inclusive
		int[]ans=new int[Math.abs(stop-start)+1];
		if (start <= stop)	for (int i=start,j=0; i<=stop; i++)		ans[j++] = i;
		else				for (int i=start,j=0; i>=stop; i--)		ans[j++] = i;
		return ans;
	}
	void transposeFaces(int w, int h) {
		// transpose 'face'
		int[][][]transpose=new int[h][w][], faces=new int[w][h][];
		for (int i=0; i<w; i++) 	for (int j=0; j<h; j++) faces[i][j]=face[i*h+j];
		for (int i=0; i<w; i++) 	for (int j=0; j<h; j++) transpose[j][i]=faces[i][j];
		for (int i=0; i<h; i++) 	for (int j=0; j<w; j++) face[i*h+j]=transpose[h-i-1][j];
	}
	public void loadBuffers(GL10 gl) {
		gl.glVertexPointer  (3, GL10.GL_FLOAT, 0, fbCoords.getBuffer());
		if (fbColor!=null) 		gl.glColorPointer  	(4, GL10.GL_FLOAT, 0, fbColor.getBuffer());
		if (fbNormal!=null) 	gl.glNormalPointer	(	GL10.GL_FLOAT, 0, fbNormal.getBuffer());
		if (fbText!=null)		gl.glTexCoordPointer(2, GL10.GL_FLOAT, 0, fbText.getBuffer());
	}
	public void drawWire(GL10 gl) {draw(gl);}
	public void drawSolidColor(GL10 gl) {
		loadBuffers(gl);
		for (int i=0,s=0; i<face.length; s+=face[i++].length) 
			gl.glDrawArrays(GL10.GL_TRIANGLE_FAN, s, face[i].length);
	}
	public void draw(GL10 gl) {
		gl.glVertexPointer(3, GL10.GL_FLOAT, 0, fbCoords.getBuffer());
		for (int i=0,s=0; i<face.length; s+=face[i++].length) 
			gl.glDrawArrays(GL10.GL_LINE_LOOP, s, face[i].length);
	}
	public void setVertex(GL10 gl) {
		gl.glVertexPointer  (3, GL10.GL_FLOAT, 0, fbCoords.getBuffer());
		gl.glTexCoordPointer(2, GL10.GL_FLOAT, 0, fbText.getBuffer());
	}
	public void drawFace(GL10 gl, int i) { // setVertex, for(getnFaces())
		gl.glDrawArrays(GL10.GL_TRIANGLE_FAN, face[0].length*i, face[i].length);
	}
	public void drawFaceDiff(GL10 gl, int i) { // when different vertex per face (pyramids)
		int s=0; for (int j=0; j<i; j++) s+=face[j].length;
		gl.glDrawArrays(GL10.GL_TRIANGLE_FAN, s, face[i].length);
	}
	public void drawSolid(GL10 gl) { // same texture in all faces
		gl.glVertexPointer  (3, GL10.GL_FLOAT, 0, fbCoords.getBuffer());
		gl.glTexCoordPointer(2, GL10.GL_FLOAT, 0, fbText.getBuffer());
		for (int i=0,s=0; i<face.length; s+=face[i++].length) 
			gl.glDrawArrays(GL10.GL_TRIANGLE_FAN, s, face[i].length);
	}
	public void drawSolidStrip(GL10 gl) { // same texture in all faces
		gl.glVertexPointer  (3, GL10.GL_FLOAT, 0, fbCoords.getBuffer());
		gl.glTexCoordPointer(2, GL10.GL_FLOAT, 0, fbText.getBuffer());
		for (int i=0,s=0; i<face.length; s+=face[i++].length) 
			gl.glDrawArrays(GL10.GL_TRIANGLE_STRIP, s, face[i].length);
	}
	public void drawFaceStrip(GL10 gl, int i) { // setVertex, for(getnFaces())
		gl.glVertexPointer  (3, GL10.GL_FLOAT, 0, fbCoords.getBuffer());
		gl.glTexCoordPointer(2, GL10.GL_FLOAT, 0, fbText.getBuffer());
		gl.glDrawArrays(GL10.GL_TRIANGLE_STRIP, face[0].length*i, face[i].length);
	}
	void scale() { // xyz -1..+1
		float max=-Float.MAX_VALUE; for (int i=0; i<xyz.length; i++) for (int j=0; j<xyz[i].length; j++) max=Math.max(max, Math.abs(xyz[i][j]));
		for (int i=0; i<xyz.length; i++) for (int j=0; j<xyz[i].length; j++) xyz[i][j]/=max;
	}
	void scaleVector(float[]v) { // v 0..+1
		float max=-Float.MAX_VALUE, min=-max, diff; 
		for (int i=0; i<v.length; i++) {
			max=Math.max(max, v[i]);	min=Math.min(min, v[i]);
		} diff=max-min;
		if (diff!=0)	for (int i=0; i<v.length; i++) v[i]=(v[i]-min)/diff;
	}
	public int getnFaces() 						{ return (face!=null) ? face.length:0;}
	public float[]getCoord(int nface, int ic) 	{ return xyz[face[nface][ic]];}
	public int getnVertexFace() 				{ return (face!=null) ? face[0].length:0;}	
	public int nVertex() { int n=0; for (int i=0; i<face.length; i++) n+=face[i].length; return n; }
	void createVertexBuffer() { // xyz[face] -> fb in face order (xyz[face[i][j])
		scale();
		fbCoords=new fbVertex(nVertex() * 3);
		for (int i=0; i<face.length; i++) 	for (int j=0; j<face[i].length; j++)	fbCoords.add(xyz[face[i][j]]); 
	}
	void createVertexBufferNoScale() { // xyz[face] -> fb in face order (xyz[face[i][j])
		fbCoords=new fbVertex(nVertex() * 3);
		for (int i=0; i<face.length; i++) 	for (int j=0; j<face[i].length; j++)	fbCoords.add(xyz[face[i][j]]); 
	}
	void createVertexBuffer(float[]vy, int nfpc) { // xyz[face] -> fb in face order (xyz[face[i][j]), nfpc=n faces per coord
		scale(); 
		fbCoords=new fbVertex(nVertex() * 3);
		for (int i=0, iv=0; i<face.length; i++, iv++)
			for (int j=0; j<face[i].length; j++) {
				for (int k=0; k<3; k++) fbCoords._add(xyz[face[i][j]][k] * vy[iv/nfpc]);
			}
		fbCoords.reset();
	}
	void createTextureBuffer() { // texture[face] -> fb in face order (texture[face[i][j])
		fbText=new fbVertex(nVertex() * 2);
		for (int i=0; i<face.length; i++) 	for (int j=0; j<face[i].length; j++)	fbText.add(texture[face[i][j]]); 
	}
	void createColorBuffer() { // color[face] -> fb in face order (color[face[i][j])
		fbColor=new fbVertex(nVertex() * 4); // rgba
		for (int i=0; i<face.length; i++) 	for (int j=0; j<face[i].length; j++)	fbColor.add(color[face[i][j]]); 
	}
	public void setColors(float[][]color) { // color[face] -> fb in face order (color[face[i][j])
		this.color=color;
		fbColor=new fbVertex(nVertex() * 4); // rgba
		for (int i=0; i<face.length; i++) 	for (int j=0; j<face[i].length; j++)	fbColor.add(color[face[i][j]]); 
	}
	void createNormalBuffer() { // each face has the same normal in every point of it (flat faces)
		fbNormal=new fbVertex(nVertex() * 3); 
		for (int i=0; i<face.length; i++) 	{
			float[]nrm=calcNormal(xyz[face[i][0]], xyz[face[i][1]], xyz[face[i][2]]);
			for (int j=0; j<face[i].length; j++)	{
				fbNormal.add(nrm);
			}
		}
	}
	float[]calcNormal(float[]p, float[]p1, float[]p2)	{ // to norm
		float xa=p1[0] - p[0], ya=p1[1] - p[1], za=p1[2] - p[2];
		float xb=p2[0] - p[0], yb=p2[1] - p[1], zb=p2[2] - p[2];
		float nx=ya * zb - za * yb, ny=za * xb - xa * zb, nz=xa * yb - ya * xb;
		return normalise(nx,ny,nz);
	}
	float[]normalise(float x, float y, float z) { // Normalise a vector
		float length = x * x + y * y + z * z;
		if (length > 0) { length = (float)Math.sqrt(length);	x /= length;	y /= length;	z /= length; } 
		else {	x = y = z = 0; }	
		return new float[]{x,y,z};
	}
	void createColorBuffer(float[]vy, int nfpc) { // color[face] -> fb in face order (color[face[i][j])
		fbColor=new fbVertex(nVertex() * 4); // rgba
		for (int i=0; i<face.length; i++) 	
			for (int j=0; j<face[i].length; j++)	fbColor.add(colIx.getColorRGBA(vy[i/nfpc])); 
	}
	void createPolyTextureBuffer() { // all poly equal nº. sides
		int ns=face[0].length;
		fbText=new fbVertex(face.length * ns * 2); // *2 = xy
		for (int i=0; i<face.length; i++) 	fbText.add(textPolygon(ns));
	}
	void createPolyTextureBufferDiff() { // all poly NOT equal nº. sides
		int ns=0;	for (int i=0; i<face.length; i++) 	ns+=face[i].length; // calc munber of sides 'ns'
		fbText=new fbVertex(ns * 2); // *2 = xy
		for (int i=0; i<face.length; i++) 	fbText.add(textPolygon(face[i].length));
	}
	float[]textPolygon(int n) {
		float[]c=new float[n*2]; 
		float pi2n=6.28318530718f/n, r=0.5f, xoffset=0.5f, yoffset=xoffset;
		for (int i=0; i<n; i++) {	float a=n-i-n/2f; c[i*2+0]=r*(float)Math.sin(a*pi2n)+xoffset; c[i*2+1]=r*(float)Math.cos(a*pi2n)+yoffset;	}
		return c;
	}
	float[]textPolygon(int n, float xoffset) {
		float[]c=new float[n*2]; 
		float pi2n=6.28318530718f/n, r=0.5f, yoffset=xoffset;
		for (int i=0; i<n; i++) {	float a=n-i-n/2f; c[i*2+0]=r*(float)Math.sin(a*pi2n)+xoffset; c[i*2+1]=r*(float)Math.cos(a*pi2n)+yoffset;	}
		return c;
	}
	//-------------------------primitive polyhedra-----------------
	public Polyhedron tetrahedron() {
		name="T";
		face=new int[][]{{0,1,2}, {0,2,3}, {0,3,1}, {1,3,2}};
		xyz=new float[][] {{1,1,1},{1,-1,-1},{-1,1,-1},{-1,-1,1}};
		createVertexBuffer();
		createPolyTextureBuffer();
		return this;
	}
	public Polyhedron tetrahedron(int col) {
		dimcolor(4);
		for (int i=0; i<4; i++) color[i]=ColorScale.color2fv(col);
		tetrahedron();
		createColorBuffer();
		createNormalBuffer();
		return this;
	}
	public Polyhedron tetrahedron(int[] col) {
		dimcolor(4);
		for (int i=0; i<4; i++) color[i]=ColorScale.color2fv(col[i]);
		tetrahedron();
		createColorBuffer();
		createNormalBuffer();
		return this;
	}
	public void octahedron() {
		name="O";
		face=new int[][]{ {0,1,2},{0,2,3},{0,3,4},{0,4,1}, {1,4,5},{1,5,2},{2,5,3},{3,5,4} };
		xyz=new float[][]{ {0,0,1.414f},{1.414f,0,0}, {0,1.414f,0},	{-1.414f,0,0},{0,-1.414f,0},{0,0,-1.414f}};
		createVertexBuffer();
		createPolyTextureBuffer();
	}
	public Polyhedron cube() {
		name="C"; float u=0.707f;
		face=new int[][]{{3,0,1,2}, {3,4,5,0}, {0,5,6,1}, {1,6,7,2}, {2,7,4,3}, {5,4,7,6} };
		xyz=new float[][]{
				{u,u,u},   {-u,u,u}, 	{-u,-u,u}, {u,-u,u}, 
				{u,-u,-u}, {u,u,-u}, 	{-u,u,-u}, {-u,-u,-u} };

		createVertexBuffer();
		createPolyTextureBuffer();
		return this;
	}
	public Polyhedron cube(int col) { // cube color
		dimcolor(8);
		for (int i=0; i<8; i++) color[i]=ColorScale.color2fv(col);
		cube();
		createColorBuffer();
		createNormalBuffer();
		return this;
	}
	public Polyhedron  cubeStrip() { // used to fix a squared texture (draw with drawSolidStrip)
		name="C"; float u=0.5f;
		face=new int[][]{{0,1,2,3}, {4,5,6,7}, {8,9,10,11}, {12,13,14,15}, {16,17,18,19}, {20,21,22,23} };
		xyz=new float[][]{ 
				{-u, -u, u}, {u, -u, u},    {-u, u, u}, {u, u, u},// FRONT
				{-u, -u, -u}, {-u, u, -u},   {u, -u, -u}, {u, u, -u},// BACK
				{-u, -u, u}, {-u, u, u},    {-u, -u, -u}, {-u, u, -u},// LEFT
				{u, -u, -u}, {u, u, -u},    { u, -u, u}, {u, u, u},// RIGHT
				{-u, u, u},  {u, u, u},     { -u, u, -u}, {u, u, -u},// TOP
				{-u, -u, u}, {-u, -u, -u},  { u, -u, u}, {u, -u, -u} };// BOTTOM
		texture = new float[][] { 
				{0, 1}, {1, 1}, {0, 0}, {1, 0},// FRONT
				{1, 1}, {1, 0}, {0, 1}, {0, 0},// BACK
				{1, 1}, {1, 0}, {0, 1}, {0, 0},// LEFT
				{1, 1}, {1, 0}, {0, 1}, {0, 0},// RIGHT
				{1, 0}, {0, 0}, {1, 1}, {0, 1},// TOP
				{0, 0}, {0, 1}, {1, 0}, {1, 1} };// BOTTOM
		createVertexBuffer();
		createTextureBuffer();
		return this;
	}

	public void icosahedron() {
		name="I";
		face=new int[][]{
				{0,1,2}, {0,2,3}, {0,3,4}, {0,4,5}, 				{0,5,1}, {1,5,7}, {1,7,6}, {1,6,2}, 
				{2,6,8}, {2,8,3}, {3,8,9}, {3,9,4}, 				{4,9,10}, {4,10,5}, {5,10,7}, {6,7,11}, 
				{6,11,8}, {7,10,11}, {8,11,9}, {9,11,10} };
		xyz=new float[][]{
				{0,0,1.176f}, {1.051f,0,0.526f}, 				{0.324f,1,0.525f}, {-0.851f,0.618f,0.526f}, 
				{-0.851f,-0.618f,0.526f}, {0.325f,-1,0.526f},	{0.851f,0.618f,-0.526f}, {0.851f,-0.618f,-0.526f},
				{-0.325f,1.f,-0.526f}, {-1.051f,0,-0.526f}, 	{-0.325f,-1.f,-0.526f}, {0,0,-1.176f} };
		createVertexBuffer();
		createPolyTextureBuffer();
	}
	public void dodecahedron() {
		name="D";
		face=new int[][]{
				{0,1,4,7,2}, {0,2,6,9,3}, {0,3,8,5,1}, 				{1,5,11,10,4}, {2,7,13,12,6}, {3,9,15,14,8}, 
				{4,10,16,13,7}, {5,8,14,17,11}, {6,12,18,15,9},		{10,11,17,19,16}, {12,13,16,19,18}, {14,15,18,19,17} };
		xyz=new float[][]{
				{0,0,1.07047f}, {0.713644f,0,0.797878f},				{-0.356822f,0.618f,0.797878f}, {-0.356822f,-0.618f,0.797878f}, 
				{0.797878f,0.618034f,0.356822f}, {0.797878f,-0.618f,0.356822f},			{-0.934172f,0.381966f,0.356822f}, {0.136294f,1.f,0.356822f}, 
				{0.136294f,-1.f,0.356822f}, {-0.934172f,-0.381966f,0.356822f}, 			{0.934172f,0.381966f,-0.356822f}, {0.934172f,-0.381966f,-0.356822f}, 
				{-0.797878f,0.618f,-0.356822f}, {-0.136294f,1.f,-0.356822f}, 			{-0.136294f,-1.f,-0.356822f}, {-0.797878f,-0.618034f,-0.356822f}, 
				{0.356822f,0.618f,-0.797878f}, {0.356822f,-0.618f,-0.797878f},			{-0.713644f,0,-0.797878f}, {0,0,-1.07047f} };
		createVertexBuffer();
		createPolyTextureBuffer();
	}
	public Polyhedron dodecahedron(int col) { // cube color
		dimcolor(20);
		for (int i=0; i<20; i++) color[i]=ColorScale.color2fv(col);
		dodecahedron();
		createColorBuffer();
		createNormalBuffer();
		return this;
	}
	public void prism(int n) { prism(n, 1, (float)Math.sin(6.283185f/n/2)); }  
	public void prism(int n, float r, float h) {         
		float theta = pi2/n;        // pie angle
		name="P" + n;

		dimxyz(n*2);
		for (int i=0; i<n; i++) addVxyz(i, 	  new float[]{r*(float)Math.cos(i*theta), r*(float)Math.sin(i*theta), h});    // vertex #"s 0...n-1 around one face
		for (int i=0; i<n; i++) addVxyz(i+n, new float[]{r*(float)Math.cos(i*theta), r*(float)Math.sin(i*theta),-h});    // vertex #"s n...2n-1 around other
		dimface(2+n);
		face[0] = sequence(n-1, 0); face[1] = sequence(n, 2*n-1);         // top, bottom
		for (int i=0; i<n; i++) face[i+2] = new int[]{i, (i+1)%n, (i+1)%n+n, i+n};  // n square sides:
		createVertexBuffer();
		createPolyTextureBuffer();
	}
	public void antiprism(int n) {         
		name="A" + n;
		float theta = pi2/n;        // pie angle
		float h = (float)Math.sqrt(1f-4f/(4+2*(float)Math.cos(theta/2)-2*(float)Math.cos(theta))); // half-height
		float r = (float)Math.sqrt(1-h*h);      // radius of face circle
		float f = (float)Math.sqrt(h*h + (float)Math.pow(r*(float)Math.cos(theta/2), 2)); 
		r=r/f;  // correction so edge midpoints (not vertices) on unit sphere
		h=h/f;

		dimxyz(n*2);
		for (int i=0; i<n; i++)     // vertex #"s 0...n-1 around one face
			addVxyz(i, new float[]{(float)Math.cos(i*theta), (float)Math.sin(i*theta), h});
		for (int i=0; i<n; i++)     // vertex #"s n...2n-1 around other
			addVxyz(i+n, new float[]{(float)Math.cos((i+.5f)*theta), (float)Math.sin((i+.5f)*theta),-h});

		dimface(n*2+2);
		face[0] = sequence(n-1, 0);      // top
		face[1] = sequence(n, 2*n-1);    // bottom
		for (int i=0; i<n; i++) {                          // 2n triangular sides:
			face[i*2+2 +0] = new int[]{i, (i+1)%n, i+n};
			face[i*2+2 +1] = new int[]{i, i+n, ((n+i-1)%n+n)};
		}
		//		xyz = adjustXYZ(1);
		createVertexBuffer();
		createPolyTextureBuffer();
	}
	public void pyramid(int n) {         
		name="Y" + n;
		float theta = pi2/n;        // pie angle
		int n1=n+1;
		dimxyz(n1);
		for (int i=0; i<n; i++)     // vertex #"s 0...n-1 around base
			addVxyz(i,new float[]{(float)Math.cos(i*theta), (float)Math.sin(i*theta), 0});
		addVxyz(n,new float[]{0, 0, -sqPhi});    // 'n' vertex -> apex (egyptian pyramid 1, sqrt(phi)=1.272, phi)

		dimface(n1);
		face[0] = sequence(n-1, 0);      // base
		for (int i=0; i<n; i++) face[i+1] = new int[]{ n, (i+1)%n, i };   // CW: n triangular sides:

		createVertexBuffer();
		createPolyTextureBufferDiff();
	}
	public Polyhedron sphere(float rad) {	return sphere(rad,30);	}
	public Polyhedron sphere(float rad, int sec) { return sphere(new float[]{0,0,0}, rad,30,new float[]{1,1,1,1}); }	
	public Polyhedron sphere(float[]center, float rad, int sec, float[]spColor) { // center=float[3]{x,y,z}
		int n = sec+1, n2=n*n, nface=sec*sec; // calc points
		dimxyz(n2); dimface(nface); dimtexture(n2); dimcolor(n2); dimnormal(n2);

		int cv=0; float v=0, u=0, du=1f/sec, dv=1f/sec;
		for (int r=0; r<=sec; r++) {
			v=r*dv; // [0,1]
			float theta1 = v * pi; // [0,PI]

			float 	xn=0, yn=1, zn=0; // n=(0,1,0)
			float 	cosRY = (float)Math.cos(theta1), sinRY = (float)Math.sin(theta1), // n.rotateZ(theta1)
					xt = -(yn*sinRY), yt = +(yn*cosRY),	zt = 0;
			xn=xt; yn=yt; zn=zt;

			for (int c = 0; c <= sec; c++) {
				u = c * du; // [0,1]
				float theta2 = u * pi2; // [0,2PI]
				float x=xn, y=yn, z=zn; // xyz=n

				cosRY = (float)Math.cos(theta2); sinRY = (float)Math.sin(theta2); // xyz.rotateY
				xt = (x*cosRY)+(z*sinRY);	zt = (x*-sinRY)+(z*cosRY);
				x=xt; z=zt;

				x*=rad; y*=rad; z*=rad; // xyz*=rad
				x+=center[0]; y+=center[1]; z+=center[2]; // xyz+=center

				xyz[cv]=new float[]{x,y,z};
				texture[cv]=new float[]{u,v};
				cv++;
			}
		}
		int cl = sec + 1, cf=0; 		// Add faces
		for (int r = 0; r < sec; r++) {
			int off = r * cl; 
			for (int c = 0; c < sec; c++) 
				face[cf++]=new int[]{off+c,	 off+c+1, off+(c + 1 + cl), off+(c + 0 + cl)};
		}
		for (int i=0; i<n2; i++) color[i]=spColor; // set colors
		createVertexBufferNoScale(); 
		createTextureBuffer();
		createColorBuffer();
		createNormalBuffer();

		return this;
	}
	public void sphereGraph(float rad, float[]vy) { // plot a 'vy' [0..1] array in a sphere, moving faces 
		int sec=(int)Math.sqrt(vy.length);
		int n = sec+1, n2=n*n, nface=sec*sec; // calc points
		dimxyz(n2); dimface(nface); 

		int cv=0; float v=0, u=0, du=1f/sec, dv=1f/sec;
		for (int r=0; r<=sec; r++) {
			v=r*dv; // [0,1]
			float theta1 = v * pi; // [0,PI]

			float 	xn=0, yn=1, zn=0; // n=(0,1,0)
			float 	cosRY = (float)Math.cos(theta1), sinRY = (float)Math.sin(theta1), // n.rotateZ(theta1)
					xt = -(yn*sinRY), yt = +(yn*cosRY),	zt = 0;
			xn=xt; yn=yt; zn=zt;

			for (int c = 0; c <= sec; c++) {
				u = c * du; // [0,1]
				float theta2 = u * pi2; // [0,2PI]
				float x=xn, y=yn, z=zn; // xyz=n

				cosRY = (float)Math.cos(theta2); sinRY = (float)Math.sin(theta2); // xyz.rotateY
				xt = (x*cosRY)+(z*sinRY);	zt = (x*-sinRY)+(z*cosRY);
				x=xt; z=zt;

				x*=rad; y*=rad; z*=rad; // xyz*=rad
				xyz[cv]=new float[]{x,y,z};
				cv++;
			}
		}
		for (int r = 0, cl = sec + 1, cf=0; r < sec; r++) 		 		// Add faces, face*=vy, color 
			for (int c = 0, off = r * cl; c < sec; c++) 
				face[cf++]=new int[]{off+c,	 off+c+1, off+(c + 1 + cl), off+(c + 0 + cl)};

		scaleVector(vy);
		createVertexBuffer(vy,1);
		createColorBuffer(vy,1);
	}
	public void sphereGraphPyramids(float rad, float[]vy) { // plot a 'vy' [0..1] array in a sphere
		int sec=(int)Math.sqrt(vy.length), nfpc=5; // 5 side pyramid
		int n = sec+1, n2=n*n, nface=sec*sec; // calc points
		dimxyz(n2+1); dimface(nface * nfpc); 


		int cv=0; float v=0, u=0, du=1f/sec, dv=1f/sec;
		xyz[cv++]=new float[]{0,0,0}; // pyramid peak

		for (int r=0; r<=sec; r++) {
			v=r*dv; // [0,1]
			float theta1 = v * pi; // [0,PI]

			float 	xn=0, yn=1, zn=0; // n=(0,1,0)
			float 	cosRY = (float)Math.cos(theta1), sinRY = (float)Math.sin(theta1), // n.rotateZ(theta1)
					xt = -(yn*sinRY), yt = +(yn*cosRY),	zt = 0;
			xn=xt; yn=yt; zn=zt;

			for (int c = 0; c <= sec; c++) {
				u = c * du; // [0,1]
				float theta2 = u * pi2; // [0,2PI]
				float x=xn, y=yn, z=zn; // xyz=n

				cosRY = (float)Math.cos(theta2); sinRY = (float)Math.sin(theta2); // xyz.rotateY
				xt = (x*cosRY)+(z*sinRY);	zt = (x*-sinRY)+(z*cosRY);
				x=xt; z=zt;

				x*=rad; y*=rad; z*=rad; // xyz*=rad
				xyz[cv++]=new float[]{x,y,z};
			}
		}
		for (int r = 0, cl = sec + 1, cf=0; r < sec; r++) 		 		// Add faces, face*=vy 
			for (int c = 1, off = r * cl; c <= sec; c++) { // from coord '1' as coor '0' is peak
				face[cf++]=new int[]{off+c,	 off+c+1, off+(c + 1 + cl), off+(c + 0 + cl)};
				face[cf++]=new int[]{0,0, off+c,	 off+c+1};
				face[cf++]=new int[]{0,	0, off+(c + 1 + cl), off+(c + 0 + cl)};
				face[cf++]=new int[]{0,	 off+c+1, off+(c + 1 + cl), 0};
				face[cf++]=new int[]{off+c,	 0, 0, off+(c + 0 + cl)};
			}

		scaleVector(vy);
		createVertexBuffer(vy,nfpc);
		createColorBuffer(vy,nfpc);
	}
	private class Number3d	{ public float x,y,z;	}
	private class Vertex3d {
		public Number3d		position = new Number3d();
		public Number3d		normal;
		public float[]fvPos() { return new float[]{position.x,position.y,position.z}; }
	}
	public void torus(float lr, float sr, int w, int h) { // w=secs in circle, lr=1, sr=.1f (1, .1, 30, 12) 
		int cv=0, cf=0, vcount = 0;

		dimxyz(4*w*h); dimtexture(4*w*h); dimface(w*h);
		float step1r = pi2 / w, step2r = pi2 / h, a1a = 0, a1b = step1r;

		for(float s=0; s<w; s++, a1a=a1b, a1b+=step1r) {
			float a2a = 0;
			float a2b = step2r;

			for(float s2=0; s2<h; s2++, a2a=a2b, a2b+=step2r) {
				Vertex3d v0 = getVertex(a1a, lr, a2a, sr), v1 = getVertex(a1b, lr, a2a, sr), 
						v2 = getVertex(a1b, lr, a2b, sr), v3 = getVertex(a1a, lr, a2b, sr);
				float ux1 = s/w, ux0 = (s+1)/w, uy0 = s2/h, uy1 = (s2+1)/h;

				xyz[cv+0]=v0.fvPos();	xyz[cv+1]=v1.fvPos();	xyz[cv+2]=v2.fvPos();	xyz[cv+3]=v3.fvPos();
				texture[cv+0]=Uv(1-ux1, uy0);	texture[cv+1]=Uv(1-ux0, uy0);	texture[cv+2]=Uv(1-ux0, uy1);	texture[cv+3]=Uv(1-ux1, uy1);
				cv+=4;
				face[cf++]=new int[]{vcount, vcount+1, vcount+2, vcount+3};

				vcount += 4;
			}
		}
		createVertexBuffer();
		createTextureBuffer();
	}
	private float[]Uv(float u, float v) {return new float[]{u,v};}
	private Vertex3d getVertex(float a1, float lr, float a2, float sr) {
		Vertex3d vertex = new Vertex3d();
		vertex.normal   = new Number3d();
		float ca1 = (float)Math.cos(a1), sa1 = (float)Math.sin(a1), ca2 = (float)Math.cos(a2), sa2 = (float)Math.sin(a2);
		float centerX = lr * ca1, centerZ = -lr * sa1;
		vertex.normal.x = ca2 * ca1;		vertex.normal.y = sa2;		vertex.normal.z = -ca2 * sa1;
		vertex.position.x = centerX + sr * vertex.normal.x;
		vertex.position.y = sr * vertex.normal.y;
		vertex.position.z = centerZ + sr * vertex.normal.z;
		return vertex;
	}
	public void tube(int n, float r, float h) {prism(n,r,h);}
	public void plane() {
		name="L";
		face=new int[][]{{3,0,1,2}};
		xyz=new float[][]{	{0.707f,0.707f,0.707f},   {-0.707f,0.707f,0.707f}, 				{-0.707f,-0.707f,0.707f}, {0.707f,-0.707f,0.707f}}; 
		createVertexBuffer();
	}
	public void polygon(int n) {
		float pi2n=6.28318530718f/n;
		dimxyz(n); dimface(1);
		for (int i=0; i<n; i++) xyz[i]=new float[]{(float)Math.sin(i*pi2n), (float)Math.cos(i*pi2n), 0};
		face[0]=sequence(0,n-1);
		createVertexBuffer();
	}
	public void polyStar(int n) {
		float pi2n=6.28318530718f/n;
		dimxyz(n); dimface(1);
		for (int i=0; i<n; i++) xyz[i]=new float[]{(float)Math.sin(2*i*pi2n), (float)Math.cos(2*i*pi2n), 0};
		face[0]=sequence(0,n-1);
		createVertexBuffer();
	}
	public void grid(int n) {
		float n1=1f/n, l=1, l2=l/2f; int nn1=n+1;
		dimxyz(nn1*2*2); dimface(nn1*2);
		for (int i=0,j=0; i<nn1; i++) {
			int i2=i*2, n2=nn1*2;
			xyz[i2+0] =new float[]{i*n1-l2, -l2, 0}; xyz[i2+1]   =new float[]{i*n1-l2,  l-l2, 0};   face[j++]=new int[]{i2,i2+1};
			xyz[i2+n2]=new float[]{l-l2, i*n1-l2, 0}; xyz[i2+n2+1]=new float[]{0-l2 , i*n1-l2, 0};  face[j++]=new int[]{i2+n2,i2+n2+1};
		}
		createVertexBuffer(); // gl.glScalef(4,7,4); gl.glRotatef(-77, 1, 0, 0); // creates a vert perspective
	}
	public Polyhedron panel(float w, float h) { 
		float p2 = w/2, l2=h/2, mp2=-p2, ml2=-l2;
		float[]txCoords	= {0,0, 1,0, 1,1, 0,1}; // {x0,y0, x1,y0, x1,y1, x0,y1}; TRIANGLE_STRIP texture coords order={01,11,00,10}
		xyz=new float[][]{{ml2,mp2,0},  {l2,mp2,0},  {l2,p2,0},  {ml2,p2,0}}; // CCW phi ratio rect;	
		face=new int[][]{{3,2,1,0}};
		fbText=new fbVertex(txCoords);
		createVertexBuffer();
		return this;
	}
	public Polyhedron panel() { 	return panel(1,1);	} // default panel
	public void multiPanel(int w, int h) { // w x h centered grid panel "xyz & texture=(00,10,11,01) TRIANGLE_FAN"
		int w1=w+1, h1=h+1, wh1=w1*h1, wh=w*h; float xd=1f/w, yd=1f/h;
		dimxyz(wh1); dimface(wh); dimtexture(wh);
		for (int i=0,c=0,cf=0; i<=w; i++) { // calc vertex(wh1) & faces(wh) & textures (wh)
			for (int j=0; j<=h; j++,c++) { 
				xyz[i*h1+j]=new float[]{xd*i-.5f, yd*j-.5f, 0}; // w1 x h1 vertex, centered 0,0 
				if (i<w & j<h) {
					float x0=xd*j, y0=yd*i, x1=x0+xd, y1=y0+yd;
					texture[cf]=new float[]{ x1,y0, x1,y1, x0,y1, x0,y0 }; // quad texture per face (same order)
					face   [cf]=new int[]  {c+h1+1, c+h1 , c, c+1}; 
					cf++; 
				}
			}
		}
		transposeFaces(w,h);
		fbText=new fbVertex(texture); // put textures in calc order
		createVertexBuffer(); // put vertex in quad face order	(gl.glRotatef(-90, 0, 0, 1);)
	}
	public void slidePanel(int w, int h) { // w x h centered grid multipanel for each panel=1 texture 
		int w1=w+1, h1=h+1, wh1=w1*h1, wh=w*h; float xd=1f/w, yd=1f/h;
		dimxyz(wh1); dimface(wh); dimtexture(wh);
		for (int i=0,c=0,cf=0; i<=w; i++) { // calc vertex(wh1) & faces(wh) & textures (wh)
			for (int j=0; j<=h; j++,c++) { 
				xyz[i*h1+j]=new float[]{xd*i-.5f, yd*j-.5f, 0}; // w1 x h1 vertex, centered 0,0 
				if (i<w & j<h) {
					float x0=xd*j, y0=yd*i, x1=x0+xd, y1=y0+yd;
					face   [cf]=new int[]  {c+h1+1, c+h1 , c, c+1}; 
					cf++; 
				}
			}
		}
		createVertexBuffer();
		createPolyTextureBuffer();
	}
	public void lemniscata() { lemniscata(70); }
	public void lemniscata(int n) {
		// plot2d([parametric, (cos(t))/(1+sin(t)^2), (cos(t)*sin(t))/(1+sin(t)^2), [t,0,2*%pi],[nticks,70]]);
		int cv=0;
		dimxyz(n); dimface(1);
		float t=0, delta=pi2/n; 
		for (int i=0; i<n; i++) {
			t=delta*i;
			float 	u=(float)Math.cos(t)/(1+sqr((float)Math.sin(t))),
					v=((float)Math.cos(t)*(float)Math.sin(t))/(1+sqr((float)Math.sin(t)));
			xyz[cv++]=new float[]{u,v,0};
		}
		face[0]=sequence(0,n-1);
		createVertexBuffer();
	}
	
	public int[]	getFace(int i) 	{ 	return face[i];  	} // {c0,c1..cn}
	public float[] 	getXYZ(int i) 	{ 	return xyz[i];		} // {x,y,z}
}