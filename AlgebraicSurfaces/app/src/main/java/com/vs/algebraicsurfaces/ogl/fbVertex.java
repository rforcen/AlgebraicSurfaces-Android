package com.vs.algebraicsurfaces.ogl;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;

public class fbVertex {
	FloatBuffer fb;		int cnt, size;
	public fbVertex(int n) 		{ alloc(n); }
	public fbVertex(float[]v)		{put(v);}
	public fbVertex(float[][]v) 	{
		int sz=0; for (float[]vv:v) sz+=vv.length; // calc sizeof v 
		alloc(sz); 
		for (float[]vv:v) add(vv);
	}
	public fbVertex() 		{ }

	public FloatBuffer alloc(int n)	 	{fb = (ByteBuffer.allocateDirect( n * 4 ).order(ByteOrder.nativeOrder())).asFloatBuffer(); cnt=0; size=n; return fb;}
	public void add(float []coords)		{for (int i=0; i<coords.length; i++) _add(coords[i]); reset(); }
	public void add(double[]coords, double scale)		{for (int i=0; i<coords.length; i++) _add((float)(coords[i]*scale)); reset(); }
	public void put(float []coords)		{fb=alloc(coords.length); fb.put(coords).position(0); }
	public void put(int i, float v) 		{fb.position(i);	fb.put(v); cnt=i;}
	public void add(float v) 				{fb.position(cnt++); fb.put(v); reset(); }
	public void add(float u, float v) 	{_add(u); _add(v); reset();}
	public void add(float x, float y, float z) 	{_add(x); _add(y); _add(z); reset();}
	public void _add(float v) 				{fb.position(cnt++); fb.put(v); } // no reset
	public void reset() 					{fb.position(0); } // must reset to position(0) before use!!
	public FloatBuffer getBuffer() 			{return fb; }
	public int getSize() 					{return size; }
	public int getCnt() 					{return cnt; }	
	public void addColor(int c)			{_add(ColorIndex.getRedf(c)); _add(ColorIndex.getGreenf(c)); _add(ColorIndex.getBluef(c));		_add(1f); }
	public byte[]getBytes()				{ // convert FloatBuffer -> byte[]
		ByteBuffer bb = ByteBuffer.allocateDirect( size * 4 );
		bb.order(ByteOrder.nativeOrder());
		bb.asCharBuffer();
		for (int i=0; i<size; i++)	bb.putFloat(fb.get(i));
		bb.position(0);
		if (bb.hasArray()) return bb.array();
		else {
			byte[]b=new byte[size*4];
			bb.get(b);
			return b;
		}
	}
	public void put(byte[]b) { // byte[]b -> fb
		ByteBuffer bb=ByteBuffer.allocateDirect( b.length ).order(ByteOrder.nativeOrder()); 
		bb.put(b); bb.position(0);
		fb=bb.asFloatBuffer();
	}
}
