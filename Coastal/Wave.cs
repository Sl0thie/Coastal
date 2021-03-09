using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Coastal
{
    /// <summary>
    /// Wave Object.
    /// </summary>
    public class Wave
    {

        private float _H;
        public float H
        {
            get { return _H; }
            set { _H = value; }
        }

        private float _T;
        public float T
        {
            get { return _T; }
            set { _T = value; }
        }

        private float _h;
        public float h
        {
            get { return _h; }
            set { _h = value; }
        }

        private float _f;
        public float f
        {
            get { return _f; }
            set { _f = value; }
        }

        private float _ω;
        public float ω
        {
            get { return _ω; }
            set { _ω = value; }
        }

        private float _k;
        public float k
        {
            get { return _k; }
            set { _k = value; }
        }

        private float _L;
        public float L
        {
            get { return _L; }
            set { _L = value; }
        }

        private float var_c;
        public float c
        {
            get { return var_c; }
            set { var_c = value; }
        }

        private float var_c_g;
        public float c_g
        {
            get { return var_c_g; }
            set { var_c_g = value; }
        }

        public Wave(float var_H, float var_T, float var_h)
        {
            //Assign properties.
            _H = var_H;
            _T = var_T;
            _h = var_h;

            //Basic wave properties.
            L = Functions.WaveLengthExact(1, T, h);
            f = Functions.AverageFrequency(T);
            ω = Functions.AngularFrequency(T);
            c = Functions.GetWaveCelerity(h, T, L);
            c_g = Functions.GetGroupVelosity(h, T, L, k, c);

        }

    }

}
