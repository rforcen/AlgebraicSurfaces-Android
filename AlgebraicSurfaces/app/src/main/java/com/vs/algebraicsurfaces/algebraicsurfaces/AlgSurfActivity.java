package com.vs.algebraicsurfaces.algebraicsurfaces;

import android.app.Activity;
import android.content.Context;
import android.os.Bundle;
import android.view.MotionEvent;
import android.widget.SeekBar;
import android.widget.TextView;

import com.vs.algebraicsurfaces.R;
import com.vs.algebraicsurfaces.ogl.OGLRenderer;

import javax.microedition.khronos.opengles.GL10;

import static com.vs.algebraicsurfaces.R.mipmap.rock;

public class AlgSurfActivity extends Activity {
    Renderer renderer;

    SeekBar seekBar; // UI items
    TextView text;

    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_alg_surf);
        initControls();
    }

    void initControls() {
        if (renderer == null) renderer = new Renderer(this, true);

        text = (TextView) findViewById(R.id.textView);
        text.setText(renderer.getName());

        seekBar = (SeekBar) findViewById(R.id.seekBar);
        seekBar.setMax(renderer.st.srfNames.length - 1);
        seekBar.setOnSeekBarChangeListener(new SeekBar.OnSeekBarChangeListener() {
            @Override
            public void onProgressChanged(SeekBar seekBar, int progress, boolean fromUser) {
                text.setText(renderer.getName(progress));
            }

            @Override
            public void onStartTrackingTouch(SeekBar seekBar) {
            }

            @Override
            public void onStopTrackingTouch(SeekBar seekBar) {
                renderer.setn(seekBar.getProgress());
            }
        });

    }

    @Override
    public boolean onTouchEvent(MotionEvent event) {
        renderer.cursorMove(event);
        renderer.refresh();
        return false;
    }

    class Renderer extends OGLRenderer { /////////////////////////// RENDERER
        Context context;
        SurfaceTest st = new SurfaceTest();
        int ns = 0, resolution = 80;

        public Renderer(Context context, boolean hasTx) {
            super(context, R.id.glSurface, hasTx);
            this.context = context;
        }

        public void setNext() {
            ns = (ns + 1) % st.srfNames.length;
            st.calcCoords(ns, resolution);
            refresh();
        }

        public void setn(int n) {
            ns = n;
            st.calcCoords(ns, resolution);
            refresh();
        }

        public String getName() {
            return st.srfNames[ns];
        }

        public String getName(int n) {
            return st.srfNames[n];
        }

        @Override
        public void drawModel(GL10 gl) {
            st.draw(gl, ns, true);
            postRotate();
        }

        @Override
        public void postCreate(GL10 gl) {
            setZoom(-3.7f);
            loadTexture(gl, rock);
            st.calcCoords(ns, resolution);
        }
    }
}
