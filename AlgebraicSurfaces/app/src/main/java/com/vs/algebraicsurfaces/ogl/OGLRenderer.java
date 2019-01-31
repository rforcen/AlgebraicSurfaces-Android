package com.vs.algebraicsurfaces.ogl;

import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.ShortBuffer;

import javax.microedition.khronos.egl.EGLConfig;
import javax.microedition.khronos.opengles.GL10;
import android.app.Activity;
import android.content.Context;
import android.graphics.Bitmap;
import android.graphics.Bitmap.CompressFormat;
import android.graphics.BitmapFactory;
import android.opengl.GLSurfaceView;
import android.opengl.GLU;
import android.opengl.GLUtils;
import android.util.DisplayMetrics;
import android.view.MotionEvent;
import android.view.ScaleGestureDetector;
import android.view.ScaleGestureDetector.OnScaleGestureListener;

/* 												minimal usage
class Renderer extends OGLRenderer { //	 		renderer
		public Renderer(Context context, boolean hasTx) {	 super(context,R.id.glSurface, hasTx); }
		@Override public void drawModel(GL10 gl) 	{}
		@Override public void postCreate(GL10 gl) 	{}
	}
 */

// 												generic OGL renderer -> implement drawModel
public abstract class OGLRenderer implements GLSurfaceView.Renderer  {
	public volatile float angleX=0, angleY=0, angleZ=0;
	public  float zoom=-3.5f, moveSpeed=.1f;
	public Context context;
	float mPreviousX=0, mPreviousY=0, TOUCH_SCALE_FACTOR=180f/1320f; // move stuff
	boolean hasTexture=false;
	private long startTime, fpsStartTime, numFrames;
	public  float fps, afps;
	public 	 GLSurfaceView 	glSurface;	
	private int 	 idResource;
	private ScaleGestureDetector mScaleDetector; // scale gesture
	private float scale=1, currentScale=1;
	private int w,h;
	private boolean screenshot=false; // screen save: call setScreenShot
	private String fnScreenShot;
	private CompressFormat ssFmt= CompressFormat.JPEG; // default ss (jpg, 90)
	private int ssquality=90;


	public OGLRenderer(Context context, int idResource, boolean hasTexture) {	
		this.context = context; this.idResource=idResource; this.hasTexture=hasTexture;
		initScaleGesture();
		setOGL();
	}
	public OGLRenderer(GLSurfaceView glSurface, boolean hasTexture) {
		this.hasTexture=hasTexture;
		glSurface.setRenderer(this);
	}
	public float 	getCurrentScale()		{ return scale*currentScale; }
	public abstract void drawModel(GL10 gl); // implement it on extended class
	public abstract void postCreate(GL10 gl); /* howto load multiple textures: 	
		gl.glGenTextures(numTextures, int[]Textures, offset=0); // generate the texture set 
    	gl.glBindTexture(GL10.GL_TEXTURE_2D, Textures[i]);
		loadTexture(gl, R.drawable.name);	*/

	private void setOGL() { // set the OGL env.
		glSurface=(GLSurfaceView)((Activity) context).findViewById(idResource); 
		glSurface.setRenderer(this);
	}

	public void initView(GL10 gl) {// Position model 
		gl.glClear(GL10.GL_COLOR_BUFFER_BIT | GL10.GL_DEPTH_BUFFER_BIT);  // Clear the screen to black
		identView(gl);
	}
	public void identView(GL10 gl) {
		gl.glMatrixMode(GL10.GL_MODELVIEW);
		gl.glLoadIdentity();
		gl.glTranslatef(0, 0, zoom);
	}
	@Override public synchronized void onDrawFrame(GL10 gl) {
		initView(gl);
		rotate(gl);
		drawModel(gl); // draw model
		saveScreenShot(gl); // triggered by setScreenShot
	}
	@Override	public void onSurfaceChanged(GL10 gl, int w, int h) {//Define the view frustum
		this.w=w; this.h=h;
		gl.glViewport(0, 0, w, h);
		gl.glMatrixMode(GL10.GL_PROJECTION);
		gl.glLoadIdentity();
		float ratio = (float) w / h;
		GLU.gluPerspective(gl, 45.0f, ratio, 1, 100f); 
	}
	@Override	public void onSurfaceCreated(GL10 gl, EGLConfig arg1) {
		startTime = System.currentTimeMillis();	fpsStartTime = startTime;	numFrames = 0;

		if (hasTexture) enableTexture(gl);
		gl.glEnable(GL10.GL_DEPTH_TEST); 
		gl.glDepthFunc(GL10.GL_LEQUAL);
		gl.glEnableClientState(GL10.GL_VERTEX_ARRAY);

		gl.glDisable(GL10.GL_DITHER);
		postCreate(gl);
	}
	public void rotate(GL10 gl) {	gl.glRotatef(angleX, 0,1,0);gl.glRotatef(angleY, 1,0,0);gl.glRotatef(angleZ, 0,0,1);}
	public boolean  cursorMove(MotionEvent event) { // updates angleX, angleY -> must render
		final DisplayMetrics displayMetrics = new DisplayMetrics();
		((Activity)context).getWindowManager().getDefaultDisplay().getMetrics(displayMetrics);
		float density = displayMetrics.density;
		TOUCH_SCALE_FACTOR=density/5;

		float 	x = event.getX(), y = event.getY();

		switch (event.getAction()) {
		case MotionEvent.ACTION_MOVE:
			float 	dx = x - mPreviousX, dy = y - mPreviousY;
			angleX += dx * TOUCH_SCALE_FACTOR;  
			angleY += dy * TOUCH_SCALE_FACTOR;  
			break;
		case MotionEvent.ACTION_DOWN:
			break;
		}
		mPreviousX = x;		mPreviousY = y;
		return true;
	}
	public boolean cursorScale(MotionEvent event) { return mScaleDetector.onTouchEvent(event);	}
	public boolean rotateScale(MotionEvent event) { if (cursorScale(event)) {	cursorMove(event);	refresh();	} return false;}
	
	public boolean lap(float secs) {
		boolean done=false;
		numFrames++;
		long fpsElapsed = System.currentTimeMillis() - fpsStartTime;

		if (fpsElapsed > secs * 1000) { // every 5 seconds
			fps = (numFrames * 1000.0F) / fpsElapsed;
			if (afps==0) { afps+=fps; } // average
			else { afps+=fps; afps/=2;}
			fpsStartTime = System.currentTimeMillis();
			numFrames = 0;
			done=true;
		}
		return done;
	}
	public void setZoom(float zoom) { this.zoom=zoom; }
	public void loadTexture(GL10 gl,  int resource) {
		Bitmap bmp = BitmapFactory.decodeResource(context.getResources(), resource);
		if (bmp!=null)	loadTexture(gl, bmp);
	}
	public void loadTexture(GL10 gl, Bitmap bmp ) {
		if (bmp==null) return;
		GLUtils.texImage2D(GL10.GL_TEXTURE_2D, 0, bmp, 0);
		gl.glTexParameterx(GL10.GL_TEXTURE_2D,  GL10.GL_TEXTURE_MIN_FILTER, GL10.GL_LINEAR);
		gl.glTexParameterx(GL10.GL_TEXTURE_2D,  GL10.GL_TEXTURE_MAG_FILTER, GL10.GL_LINEAR);
		bmp.recycle();
	}
	public Bitmap getBitmapFromAsset(String strName) {
		Bitmap bitmap = null;
		try { bitmap = BitmapFactory.decodeStream(context.getAssets().open(strName)); } catch (IOException e) {}
		return bitmap;
	}
	public void loadTextureAsset(GL10 gl, String fileAsset ) {	loadTexture(gl, getBitmapFromAsset(fileAsset));	}
	public void light(GL10 gl) {// create light
		//Define the lighting
		float lightAmbient[] = new float[] { 0.2f, 0.2f, 0.2f, 1 };
		float lightDiffuse[] = new float[] { 1, 1, 1, 1 };
		float[] lightPos0 = new float[] { 1, 1, 1, 1 };
		float[] lightPos1 = new float[] { -5, -5, -5, 0};
		gl.glEnable(GL10.GL_LIGHTING);
		gl.glEnable(GL10.GL_LIGHT0);
		gl.glLightfv(GL10.GL_LIGHT0, GL10.GL_AMBIENT, lightAmbient, 0);
		gl.glLightfv(GL10.GL_LIGHT0, GL10.GL_DIFFUSE, lightDiffuse, 0);
		gl.glLightfv(GL10.GL_LIGHT0, GL10.GL_POSITION, lightPos0, 0);

		gl.glEnable(GL10.GL_LIGHT1);
		gl.glLightfv(GL10.GL_LIGHT1, GL10.GL_AMBIENT, lightAmbient, 0);
		gl.glLightfv(GL10.GL_LIGHT1, GL10.GL_DIFFUSE, lightDiffuse, 0);
		gl.glLightfv(GL10.GL_LIGHT1, GL10.GL_POSITION, lightPos1, 0);
	}
	public void material(GL10 gl) {
		float matAmbient[] = new float[] { 1, 1, 1, 1 };
		float matDiffuse[] = new float[] { 1, 1, 1, 1 };
		gl.glMaterialfv(GL10.GL_FRONT_AND_BACK, GL10.GL_AMBIENT,matAmbient, 0);
		gl.glMaterialfv(GL10.GL_FRONT_AND_BACK, GL10.GL_DIFFUSE,matDiffuse, 0);
	}
	public void enableTexture(GL10 gl) {
		light(gl); material(gl); 
		gl.glEnableClientState(GL10.GL_TEXTURE_COORD_ARRAY);
		gl.glEnable(GL10.GL_TEXTURE_2D);
	}
	public void setTexture(GL10 gl, int nt) {
		gl.glBindTexture(GL10.GL_TEXTURE_2D, nt);
	}
	public void startAnimation() 	{ glSurface.setRenderMode(GLSurfaceView.RENDERMODE_CONTINUOUSLY); }
	public void stopAnimation() 	{ glSurface.setRenderMode(GLSurfaceView.RENDERMODE_WHEN_DIRTY); }
	public void switchAnimation() 	{ glSurface.setRenderMode((glSurface.getRenderMode()==GLSurfaceView.RENDERMODE_CONTINUOUSLY) ? GLSurfaceView.RENDERMODE_WHEN_DIRTY:GLSurfaceView.RENDERMODE_CONTINUOUSLY); }
	public void refresh()			{ glSurface.requestRender(); }
	public void postRotate() 		{ if (glSurface.getRenderMode()==GLSurfaceView.RENDERMODE_CONTINUOUSLY)  {angleX+=moveSpeed;  angleY+=moveSpeed; angleZ+=moveSpeed;} }
	public void setMoveSpeed(float moveSpeed) 		{	this.moveSpeed=moveSpeed;}
	public void enableColorArray(GL10 gl) 			{	gl.glEnableClientState(GL10.GL_COLOR_ARRAY); }
	public void enableColorArrayMaterial(GL10 gl)	{	gl.glEnableClientState(GL10.GL_COLOR_ARRAY); gl.glEnable(GL10.GL_COLOR_MATERIAL); }
	public void enableNormalArray(GL10 gl)			{	gl.glEnableClientState(GL10.GL_NORMAL_ARRAY); }
	public void setScale(GL10 gl, float scale)		{	gl.glScalef(scale, scale, scale);	}
	public void enableTransparency(GL10 gl) 		{	gl.glEnable(GL10.GL_BLEND); }
	public void setColor(GL10 gl, int color) 						{	gl.glColor4f(ColorScale.getRedf(color), ColorScale.getGreenf(color), ColorScale.getBluef(color), 1);	}
	public void setColor(GL10 gl, int color, float alpha) 		{	gl.glColor4f(ColorScale.getRedf(color), ColorScale.getGreenf(color), ColorScale.getBluef(color), alpha);	}
	public void disableTextures(GL10 gl) 			{ gl.glDisable(GL10.GL_TEXTURE_2D); }	
	public void enableTextures(GL10 gl) 			{ gl.glEnable(GL10.GL_TEXTURE_2D); }
	public void resetScale()						{ scale=currentScale=1; refresh(); }

	private void initScaleGesture() {
		mScaleDetector = new ScaleGestureDetector(context, new OnScaleGestureListener() {
			@Override public void 		onScaleEnd(ScaleGestureDetector detector) 	{ fixScale(); }
			@Override public boolean 	onScaleBegin(ScaleGestureDetector detector) { return true;}
			@Override public boolean 	onScale(ScaleGestureDetector detector) 		{
				scaleModel(detector.getScaleFactor());
				return false;
			}
		});
	}
	private void 	scaleModel(float scale)	{ this.scale=scale; refresh(); }
	private void 	fixScale() 				{ currentScale*=scale; scale=1; }

	public void setScreenShot(String fnScreenShot) {this.fnScreenShot=fnScreenShot; screenshot=true; refresh();}
	public void saveScreenShot(GL10 gl) { // take a SS in file defined by 'setScreenShot'
		if(screenshot) {                     
			int screenshotSize = w * h;
			ByteBuffer bb = ByteBuffer.allocateDirect(screenshotSize * 4).order(ByteOrder.nativeOrder());
			gl.glReadPixels(0, 0, w, h, GL10.GL_RGBA, GL10.GL_UNSIGNED_BYTE, bb);
			int pixelsBuffer[] = new int[screenshotSize];
			bb.asIntBuffer().get(pixelsBuffer);	bb = null;
			Bitmap bitmap = Bitmap.createBitmap(w, h, Bitmap.Config.RGB_565);
			bitmap.setPixels(pixelsBuffer, screenshotSize-w, -w, 0, 0, w, h);
			pixelsBuffer = null;

			short sBuffer[] = new short[screenshotSize];
			ShortBuffer sb = ShortBuffer.wrap(sBuffer);
			bitmap.copyPixelsToBuffer(sb);
			for (int i = 0; i < screenshotSize; ++i) {  //Making created bitmap (from OpenGL points) compatible with Android bitmap                
				short v = sBuffer[i];
				sBuffer[i] = (short) (((v&0x1f) << 11) | (v&0x7e0) | ((v&0xf800) >> 11));
			}
			sb.rewind();
			bitmap.copyPixelsFromBuffer(sb);

			// save bitmap to screenshot file
			try {
				Bitmap.createBitmap(bitmap, 0, 0, w, h).compress(ssFmt, ssquality, new FileOutputStream(fnScreenShot));
			} catch (Exception e) {	}
			screenshot = false;
		}
	}
	public Bitmap getScreenShot(GL10 gl) { // take a SS in file defined by 'setScreenShot'
		int screenshotSize = w * h;
		ByteBuffer bb = ByteBuffer.allocateDirect(screenshotSize * 4).order(ByteOrder.nativeOrder());
		gl.glReadPixels(0, 0, w, h, GL10.GL_RGBA, GL10.GL_UNSIGNED_BYTE, bb);
		int pixelsBuffer[] = new int[screenshotSize];
		bb.asIntBuffer().get(pixelsBuffer);	bb = null;
		Bitmap bitmap = Bitmap.createBitmap(w, h, Bitmap.Config.RGB_565);
		bitmap.setPixels(pixelsBuffer, screenshotSize-w, -w, 0, 0, w, h);
		pixelsBuffer = null;

		short sBuffer[] = new short[screenshotSize];
		ShortBuffer sb = ShortBuffer.wrap(sBuffer);
		bitmap.copyPixelsToBuffer(sb);
		for (int i = 0; i < screenshotSize; ++i) {  //Making created bitmap (from OpenGL points) compatible with Android bitmap                
			short v = sBuffer[i];
			sBuffer[i] = (short) (((v&0x1f) << 11) | (v&0x7e0) | ((v&0xf800) >> 11));
		}
		sb.rewind();
		bitmap.copyPixelsFromBuffer(sb);

		return bitmap;
	}
	public void setScreenShotType(CompressFormat fmt, int quality) { this.ssFmt=fmt; this.ssquality=quality;	}
	public void sceneInit(GL10 gl, int color) { // works nice for golden solid colors (requires normals)
		float lmodel_ambient[] = {0, 0, 0, 0};
		float lmodel_twoside[] = {GL10.GL_FALSE};
		float light0_ambient[] = {0.1f, 0.1f, 0.1f, 1.0f};
		float light0_diffuse[] = {1.0f, 1.0f, 1.0f, 0.0f};
		float light0_position[] = {0.8660254f, 0.5f, 1, 0};
		float light0_specular[] = {1, 1, 1, 0};
		float bevel_mat_ambient[] = {0, 0, 0, 1};
		float bevel_mat_shininess[] = {40};
		float bevel_mat_specular[] = {1, 1, 1, 0};
		float bevel_mat_diffuse[] = {1, 0, 0, 0};

//		gl.glEnable(GL10.GL_CULL_FACE);
//		gl.glCullFace(GL10.GL_BACK);
//		gl.glEnable(GL10.GL_DEPTH_TEST);
//		gl.glClearDepthf(1);

		gl.glClearColor(ColorScale.getRedf(color), ColorScale.getGreenf(color), ColorScale.getBluef(color), 1);

		gl.glLightfv(GL10.GL_LIGHT0, GL10.GL_AMBIENT, light0_ambient, 0);
		gl.glLightfv(GL10.GL_LIGHT0, GL10.GL_DIFFUSE, light0_diffuse, 0);
		gl.glLightfv(GL10.GL_LIGHT0, GL10.GL_SPECULAR, light0_specular, 0);
		gl.glLightfv(GL10.GL_LIGHT0, GL10.GL_POSITION, light0_position, 0);
		gl.glEnable(GL10.GL_LIGHT0);

		gl.glLightModelfv(GL10.GL_LIGHT_MODEL_TWO_SIDE, lmodel_twoside, 0);
		gl.glLightModelfv(GL10.GL_LIGHT_MODEL_AMBIENT, lmodel_ambient, 0);
		gl.glEnable(GL10.GL_LIGHTING);

		gl.glMaterialfv(GL10.GL_FRONT, GL10.GL_AMBIENT, bevel_mat_ambient,0);
		gl.glMaterialfv(GL10.GL_FRONT, GL10.GL_SHININESS, bevel_mat_shininess,0);
		gl.glMaterialfv(GL10.GL_FRONT, GL10.GL_SPECULAR, bevel_mat_specular,0);
		gl.glMaterialfv(GL10.GL_FRONT, GL10.GL_DIFFUSE, bevel_mat_diffuse,0);

		gl.glEnable(GL10.GL_COLOR_MATERIAL);
//		gl.glShadeModel(GL10.GL_SMOOTH);
	}

}

