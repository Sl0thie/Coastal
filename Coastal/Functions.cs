using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;


namespace Coastal
{
    /// <summary>
    /// Coastal Functions
    /// </summary>
    /// <remarks>
    /// Naming convention. Due to the linear nature of programming some convention
    /// for superscripts and subscripts has to used. Properties all start with standard
    /// script, if a _ is encountered then the script moves to subscript. If __ is
    /// encountered then the script moves to superscript.
    /// </remarks>
    public static class Functions
    {

        #region Constants

        /// <summary>
        /// Gravity, or gravitation, is a natural phenomenon by which all things with 
        /// mass are brought toward (or gravitate toward) one another.
        /// 
        /// [Symbol: g ]
        /// [Type  : Constant ]
        /// [Unit  : m/s² ]
        /// </summary>
        /// <value>
        /// 
        /// </value>
        /// <remarks>
        /// This is a standardised value of standard gravity. While gravity varies 
        /// depending on the location at which it is taken this value is used.
        /// The variance range is currently unknown.
        /// </remarks>       
        internal const float g = 9.80665f;

        /// <summary>
        /// π is a mathematical constant, the ratio of a circle's circumference to its diameter.
        /// 
        /// [Symbol: π, pi ]
        /// [Type  : Constant ]
        /// [Unit  : Unitless ]
        /// </summary>
        /// <value>
        /// 
        /// </value>
        /// <remarks>
        /// 
        /// </remarks>       
        internal const float π = 3.1415927410125732421875f;

        /// <summary>
        /// Density of Sea Water.
        /// 
        /// [Symbol: ρ ]
        /// [Type  : density ]
        /// [Unit  : kg/m3 ]
        /// </summary>
        /// <value>
        /// 
        /// </value>
        /// <remarks>
        /// 
        /// </remarks>       
        internal const float ρ = 1035f;

        /// <summary>
        /// Density of Air.
        /// 
        /// [Symbol: ρ_air ]
        /// [Type  : density ]
        /// [Unit  : kg/m3 ]
        /// </summary>
        /// <value>
        /// 
        /// </value>
        /// <remarks>
        /// 
        /// </remarks>       
        internal const float ρ_air = 1.3f;

        #endregion

        #region Functions

        #region New




        #endregion
        
        #region Basic Wave Properties

        /// <summary>
        /// f, average frequency is number of periods per second.
        /// 
        /// [Symbol: f ]
        /// [Type  : count ]
        /// [Unit  : Hz ]
        /// </summary>
        /// <value>
        /// 
        /// </value>
        /// <remarks>
        /// 
        /// </remarks>       
        /// <param name="T"> Wave period.</param>      
        public static float AverageFrequency(float T)
        {
            Check_T(T);

            return 1f / T;
        }

        /// <summary>
        /// ω, angular frequency is radians per second.
        /// 
        /// [Symbol: ω, omega ]
        /// [Type  : angle ]
        /// [Unit  : rad ]
        /// </summary>
        /// <value>
        /// 
        /// </value>
        /// <remarks>
        /// 
        /// </remarks>   
        /// <param name="T"> Wave period.</param>
        /// <returns></returns>
        public static float AngularFrequency(float T)
        {
            Check_T(T);

            return (π * 2f) / T;

        }

        /// <summary>
        /// k, wave number is the spatial frequency of the wave in radians per metre.
        /// 
        /// [Symbol: k ]
        /// [Type  : count ]
        /// [Unit  : unitless ]
        /// </summary>
        /// <value>
        /// 
        /// </value>
        /// <remarks>
        /// 
        /// </remarks>   
        /// <param name="L">Wave Length</param>
        /// <returns></returns>
        public static float WaveNumber(float L)
        {
            Check_L(L);

            return (π * 2) / L;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="ref_L"></param>
        /// <param name="T"></param>
        /// <param name="h"></param>
        /// <returns></returns>
        public static float WaveLengthExact(float ref_L, float T, float h)// <<<<<<<<<<<<<<<<
        {
            //? Replace the ref_L with a static of the last used to see if it is faster.


            Check_L(ref_L);
            Check_T(T);
            Check_h(h);


            float lastL = ref_L;
            float nextL = float.MaxValue;
            float diff = 1;

            while (diff > 0.00001)
            {
                nextL = ((g * T * T) / (2f * π)) * Tanh((2f * π) / (h));
                diff = Math.Abs(nextL - lastL);
                //Debug.WriteLine("diff:" + diff);
                lastL = nextL;
                //Debug.WriteLine("nextL:" + nextL);
            }

            return nextL;

        }
        
        #endregion

        /// <summary>
        /// φ, phase is the position of a point in time (an instant) on a waveform cycle. 
        /// 
        /// [Symbol: φ, psi ]
        /// [Type  :  ]
        /// [Unit  :  ]
        /// </summary>
        /// <value>
        /// 
        /// </value>
        /// <remarks>
        /// 
        /// </remarks>   
        /// <param name="H"> Wave Height</param>
        /// <param name="L"> Wave Length</param>
        /// <param name="ω"> Anglular Velocity</param>
        /// <param name="k"> Wave Number</param>
        /// <param name="z"> Depth co-ordinate</param>
        /// <param name="h"> Depth to the seafloor</param>
        /// <param name="x"> Distance co-ordinate</param>
        /// <param name="t"> Time</param>
        /// <returns></returns>
        public static float SineWaveVelocityPotential(float H, float L, float ω, float k, float z, float h, float x, float t)
        {
            Check_H(H);
            Check_L(L);
            Check_ω(ω);
            Check_k(k);
            Check_z(z);
            Check_h(h);
            Check_x(x);
            Check_t(t);

            return (-g / ω) * (H / 2) * ((float)Math.Cosh(k * (z + h)) / (float)Math.Cosh(k * h)) * ((float)Math.Sin((k * x) - (ω * t)));
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="H"></param>
        /// <param name="k"></param>
        /// <param name="x"></param>
        /// <param name="ω"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        public static float SinusoidalWaveShape(float H, float k, float x, float ω, float t)
        {
            Check_H(H);
            Check_k(k);
            Check_x(x);
            Check_ω(ω);
            Check_t(t);

            return (H / 2) * (float)Math.Cos((k * x) - (ω * t));

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="k"></param>
        /// <param name="z"></param>
        /// <param name="h"></param>
        /// <returns></returns>
        public static float PressureCoefficient(float k, float z, float h)
        {
            Check_k(k);
            Check_z(z);
            Check_h(h);

            return (Cosh(k * (z + h))) / (Cosh(k * h));

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="H_avg"></param>
        /// <param name="H_s"></param>
        /// <returns></returns>
        public static float EstimateSignificantWaveHeight(float H_avg, float H_s)
        {
            Check_H_avg(H_avg);
            Check_H_s(H_s);

            float H_rms = H_avg / H_s;
            return H_s * H_rms;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="H_rms"></param>
        /// <param name="n"></param>
        /// <param name="N"></param>
        /// <returns></returns>
        public static float HeightProbibility(float H_rms, float n, float N)
        {
            Check_H_rms(H_rms);
            Check_n(n);
            Check_N(N);

            return H_rms * Sqrt(-Ln(n / N));

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="z"></param>
        /// <returns></returns>
        public static float HydrostaticPressure(float z)
        {
            Check_z(z);

            return -ρ * g * z;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="η"></param>
        /// <param name="K_p"></param>
        /// <returns></returns>
        public static float DynamicPressure(float η, float K_p)
        {
            Check_η(η);
            Check_K_p(K_p);

            return ρ * g * η * K_p;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static float FindDepthFromHydrostaticPressure(float ps)
        {
            Check_ps(ps);

            return ps / -ρ / g;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="p_max"></param>
        /// <param name="z"></param>
        /// <param name="K_p"></param>
        /// <returns></returns>
        public static float FindWaveη(float p_max, float z, float K_p)
        {
            Check_P_max(p_max);
            Check_z(z);
            Check_K_p(K_p);

            return (p_max + (ρ * g * z)) / (ρ * g * K_p);

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="k"></param>
        /// <param name="h"></param>
        /// <param name="z"></param>
        /// <param name="Ap"></param>
        /// <returns></returns>
        public static float FindWaveHeightFromAp(float k, float h, float z, float Ap)
        {
            Check_k(k);
            Check_h(h);
            Check_z(z);
            Check_Ap(Ap);

            return (2 * Ap * Sinh(k * h)) / (Cosh(k * (z + h)));
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="H"></param>
        /// <returns></returns>
        public static float EnergyDensity(float H)
        {
            Check_H(H);

            return (ρ * g * Pow(H, 2)) / 8;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="ω"></param>
        /// <param name="H"></param>
        /// <param name="k"></param>
        /// <param name="z"></param>
        /// <param name="h"></param>
        /// <returns></returns>
        public static float[] VelocityAmplitudes(float ω, float H, float k, float z, float h)
        {
            Check_ω(ω);
            Check_h(h);
            Check_k(k);
            Check_z(z);
            Check_h(h);

            float u_hat = ω * (H / 2) * ((Cosh(k * (z + h))) / (Sinh(k * h)));
            float w_hat = ω * (H / 2) * ((Sinh(k * (z + h))) / (Sinh(k * h)));
            return new float[2] { u_hat, w_hat };
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="ω"></param>
        /// <param name="H"></param>
        /// <param name="k"></param>
        /// <param name="z"></param>
        /// <param name="h"></param>
        /// <param name="x"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        public static float[] ParticleVelocites(float ω, float H, float k, float z, float h, float x, float t)
        {
            Check_ω(ω);
            Check_H(H);
            Check_h(h);
            Check_k(k);
            Check_z(z);
            Check_h(h);
            Check_x(x);
            Check_t(t);

            float u = ω * (H / 2) * ((Cosh(k * (z + h))) / (Sinh(k * h))) * (Cos((k * x) - (ω * t)));
            float w = ω * (H / 2) * ((Sinh(k * (z + h))) / (Sinh(k * h))) * (Sin((k * x) - (ω * t)));          
            return new float[2] { u, w };
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="ω"></param>
        /// <param name="H"></param>
        /// <param name="k"></param>
        /// <param name="z"></param>
        /// <param name="h"></param>
        /// <param name="x"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        public static float[] ParticleDisplacementAmplitudes(float ω, float H, float k, float z, float h, float x, float t)
        {
            Check_ω(ω);
            Check_H(H);
            Check_h(h);
            Check_k(k);
            Check_z(z);
            Check_h(h);
            Check_x(x);
            Check_t(t);

            float α_disp = (H / 2) * ((Cosh(k * (z + h))) / (Sinh(k * h))) * (Cos((k * x) - (ω * t)));
            float β_disp = (H / 2) * ((Sinh(k * (z + h))) / (Sinh(k * h))) * (Sin((k * x) - (ω * t)));
            return new float[2] { α_disp, β_disp };           
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="isDeepWater"></param>
        /// <param name="c"></param>
        /// <param name="k"></param>
        /// <param name="h"></param>
        /// <returns></returns>
        public static float GetGroupVelosity(float h, float T, float L, float k, float c)
        {
            Check_c(c);
            Check_k(k);
            Check_h(h);

            if (h / L < 0.5)
            {
                return c / 2f;
            }
            else if (h / L > 0.05)
            {
                return Sqrt(g * h);
            }
            else
            {
                return 0.5f * (1f + ((2f * k * h) / (Sinh(2f * k * h)))) * (c);
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="isDeepWater"></param>
        /// <param name="T"></param>
        /// <param name="L"></param>
        /// <returns></returns>
        public static float GetWaveCelerity(float h, float T, float L)
        {
            Check_h(h);
            Check_T(T);
            Check_L(L);

            if (h/L < 0.5)
            {
                return (g * T) / (2f * π);
            }
            else if (h / L > 0.05)
            {
                return Sqrt(g * h);
            }
            else
            {              
                return L / T;
            }           
        }

        /// <summary>
        /// Gets the wave length. (L)
        /// </summary>
        /// <param name="isDeepWater"></param>
        /// <param name="T">Wave Period</param>
        /// <param name="h">Depth</param>
        /// <returns></returns>
        public static float GetWaveLength(bool isDeepWater, float T, float h)
        {

            if (isDeepWater)
            {
                Check_T(T);

                return ((g * Pow(T, 2)) / (2f * π));
            }
            else
            {
                Check_T(T);
                Check_h(h);

                return ((g * T * T) / (2f * π)) * Pow(Tanh(Pow((((2f * π) / T) * Sqrt(h / g)), (3f / 2f))), (2f / 3f));
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="xtv"></param>
        /// <param name="H_rms"></param>
        /// <returns></returns>
        public static float ProbabilityOfAWaveHeight(float xtv, float H_rms)
        {
            Check_xtv(xtv);
            Check_H_rms(H_rms);

            return Exp(-(xtv / H_rms) * -(xtv / H_rms));

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="H_rms"></param>
        /// <param name="n"></param>
        /// <param name="N"></param>
        /// <returns></returns>
        public static float WaveHeightExceedance(float H_rms, float n, float N)
        {
            Check_H_rms(H_rms);
            Check_N(N);
            Check_N(n);

            return H_rms * Sqrt((-Log(n / N)));
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="H_rms"></param>
        /// <param name="N"></param>
        /// <returns></returns>
        public static float FindH_max(float H_rms, float N)
        {

            Check_H_rms(H_rms);
            Check_N(N);

            return H_rms * (Log(N));

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="H_rms"></param>
        /// <returns></returns>
        public static float FindH_sig(float H_rms)
        {
            Check_H_rms(H_rms);

            return Sqrt(2f) * H_rms;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="H_rms"></param>
        /// <returns></returns>
        public static float find_H_bar(float H_rms)
        {
            Check_H_rms(H_rms);

            return 0.886f * H_rms;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="H_rms"></param>
        /// <param name="xtv"></param>
        /// <returns></returns>
        public static float FindH_per(float H_rms, float xtv)
        {
            Check_H_rms(H_rms);
            Check_xtv(xtv);

            return H_rms * Sqrt(-Log(xtv / 100));

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="H_s"></param>
        /// <returns></returns>
        public static float FindH_rms(float H_s)
        {
            Check_H_s(H_s);

            return H_s / 1.42f;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="T"></param>
        /// <param name="σ"></param>
        /// <param name="ξ"></param>
        /// <param name="μ"></param>
        /// <returns></returns>
        public static float GeneralisedValueDistribution(float T, float σ, float ξ, float μ)
        {
            Check_T(T);
            Check_σ(σ);
            Check_ξ(ξ);
            Check_μ(μ);

            float y_T = -Log(1f - (1f / T));

            if (ξ == 0)
            {
                return μ - (σ * Log(y_T));
            }
            else
            {
                return μ - ((σ / ξ) * (1f - Pow(y_T, -ξ)));
            }



        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="NeHu"></param>
        /// <param name="N"></param>
        /// <param name="σ"></param>
        /// <param name="m"></param>
        /// <param name="T"></param>
        /// <param name="u"></param>
        /// <param name="ξ"></param>
        /// <returns></returns>
        public static float GeneralisedParetoDistribution(float NeHu, float N, float σ, float m, float T, float u, float ξ)
        {

            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            float ζ_u = NeHu / N;

            if (ξ == 0)
            {
                return u + (σ * Log(m * T * ζ_u));
            }
            else
            {
                return u + ((σ / ξ) * (Pow((m * T * ζ_u), ξ) - 1));
            }
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="t_d"></param>
        /// <param name="F"></param>
        /// <param name="U_10"></param>
        /// <returns></returns>
        public static float JONSWAP(float t_d, float F, float U_10)
        {
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            // Requires F,t,U_10
            // upper limits H*mo 0.243, T*p 8.13
            float tstar = (g * t_d) / (U_10);
            float Fstar = (g * F) / (U_10 * U_10);
            float Fstar_eff = Pow((tstar / 68.8f), (3f / 2f));
            //fprintf('Data = //f\n', Fstar_eff);
            float F_eff = 0;

            if (Fstar > Fstar_eff)
            {
                F_eff = Fstar_eff;
            }
            else
            {
                F_eff = Fstar;
            }


            float Hstar_mo = 0.0016f * Sqrt(F_eff);
            float H_mo = (Hstar_mo * (U_10 * U_10)) / (g);
            float Tstar_p = 0.286f * Pow(F_eff, (1f / 3f));
            //T_p = (Tstar_p * U_10) / g;
            //T_p = T_p;
            //H_mo = H_mo;

            return 0;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c"></param>
        /// <param name="ref_c"></param>
        /// <param name="θ"></param>
        /// <returns></returns>
        public static float snells_law_for_THETA_1(float c, float ref_c, float θ)
        {

            Check_c(c);
            Check_c(ref_c);
            Check_θ(θ);

            //THETA_1 = 1/asind(c/ref_c)/sind(THETA);
            return Asind(1 / ((c / ref_c) / Sind(θ)));
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c"></param>
        /// <param name="ref_c"></param>
        /// <param name="ref_θ"></param>
        /// <returns></returns>
        public static float SnellsLaw(float c, float ref_c, float ref_θ)
        {

            Check_c(c);
            Check_c(ref_c);
            Check_θ(ref_θ);

            return Asind(Sind(ref_θ) * (c / ref_c));

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="θ"></param>
        /// <param name="ref_θ"></param>
        /// <returns></returns>
        public static float RefractionCoefficent(float θ, float ref_θ)
        {
            Check_θ(θ);
            Check_θ(ref_θ);

            return Sqrt(Cosd(ref_θ) / Cosd(θ));

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="c_g"></param>
        /// <param name="ref_c_g"></param>
        /// <returns></returns>
        public static float ShoalingCoefficent(float c_g, float ref_c_g)
        {

            Check_c(c_g);
            Check_c(ref_c_g);


            return Sqrt(ref_c_g / c_g);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="ref_H"></param>
        /// <param name="K_r"></param>
        /// <param name="K_s"></param>
        /// <returns></returns>
        public static float FindH(float ref_H, float K_r, float K_s)
        {

            Check_H(ref_H);
            Check_K_r(K_r);
            Check_K_s(K_s);

            return ref_H * K_r * K_s;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="H_d"></param>
        /// <param name="H_i"></param>
        /// <returns></returns>
        public static float DiffrationCoefficient(float H_d, float H_i)
        {
            Check_H(H_d);
            Check_H(H_i);


            return H_d / H_i;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="H_r"></param>
        /// <param name="H_i"></param>
        /// <returns></returns>
        public static float RefractionCoefficient(float H_r, float H_i)
        {
            Check_H(H_r);
            Check_H(H_i);

            return H_r / H_i;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="tanβ"></param>
        /// <param name="H_o"></param>
        /// <param name="Lo"></param>
        /// <returns></returns>
        public static float SurfSimilarityParameter(float tanβ, float H_o, float Lo)
        {
            Check_tanβ(tanβ);
            Check_H(H_o);
            Check_L(Lo);


            float ξ_o = tanβ / Sqrt(H_o / Lo);
            if (ξ_o >= 4)
            {
                Debug.WriteLine("Surfzone is considered reflective");//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            }
            else
            {
                Debug.WriteLine("Surfzone is considered dissipative");
            }
            return ξ_o;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="H_b"></param>
        /// <param name="β"></param>
        /// <param name="θ_b"></param>
        /// <returns></returns>
        public static float LongshoreCurrent(float H_b, float β, float θ_b)
        {

            Check_H(H_b);
            Check_β(β);
            Check_θ(θ_b);


            float some_const = 1f;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            return some_const * Sqrt(g * H_b) * Tan(β) * Sin(2f * θ_b);

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="θ_o"></param>
        /// <returns></returns>
        public static float NielsensExplicitApproximations_K_r(float θ_o)
        {
            Check_θ(θ_o);

            return Sqrt(Cosd(θ_o));

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="H_rms"></param>
        /// <returns></returns>
        public static float NielsensExplicitApproximations_K_s(float H_rms)
        {
            Check_H_rms(H_rms);

            //return 1 / (nthroot((4 * k * h), 4)); <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            return 0;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="H_o"></param>
        /// <param name="θ_o"></param>
        /// <param name="γ_b"></param>
        /// <param name="k"></param>
        /// <returns></returns>
        public static float Find_H_b(float H_o, float θ_o, float γ_b, float k)
        {
            Check_H(H_o);
            Check_θ(θ_o);
            Check_γ(γ_b);
            Check_k(k);


            //H_b = H_o * Math.Sqrt(Math.Cosd(THETA_o)) * (1/(nthroot((4 * k * (H_o / GAMMA_b)),4)));
            return Pow(((Pow(H_o, 4f) * Cosd(θ_o) * Cosd(θ_o)) / (4 * (k / γ_b))), (1 / 5));

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="tanβ"></param>
        /// <param name="H_rms"></param>
        /// <param name="Lo"></param>
        /// <returns></returns>
        public static float FindVerticalLengthScale(float tanβ, float H_rms, float Lo)
        {
            Check_tanβ(tanβ);
            Check_H_rms(H_rms);
            Check_L(Lo);

            if (tanβ > 0.1f)
            {
                return 0.6f * tanβ * Sqrt(H_rms * Lo);
            }
            else
            {
                return 0.06f * Sqrt(H_rms * Lo);
            }

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="L_r"></param>
        /// <param name="n"></param>
        /// <param name="N"></param>
        /// <param name="z_100"></param>
        /// <returns></returns>
        public static float NaturalWavesRayleighDistributionNumber(float L_r, float n, float N, float z_100)
        {
            Check_L(L_r);
            Check_n(n);
            Check_N(N);
            Check_z_100(z_100);


            return L_r * Sqrt(-Log(n / N)) + z_100;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="R_n"></param>
        /// <param name="L_r"></param>
        /// <returns></returns>
        public static float NaturalWavesRayleighDistributionProbability(float R_n, float L_r)
        {
            Check_R_n(R_n);
            Check_L(L_r);

            return Exp(-(R_n / L_r) * -(R_n / L_r));
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="C_10"></param>
        /// <param name="U_10"></param>
        /// <param name="ψ"></param>
        /// <returns></returns>
        public static float Find_τ_w(float C_10, float U_10, float ψ)
        {
            Check_C_10(C_10);
            Check_U_10(U_10);
            Check_ψ(ψ);


            return 0.5f * ρ_air * C_10 * (U_10 * Cosd(ψ)) * (U_10 * Cosd(ψ));

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="αν"></param>
        /// <param name="τ_w"></param>
        /// <param name="refh"></param>
        /// <param name="h"></param>
        /// <param name="W"></param>
        /// <returns></returns>
        public static float WindSetup(float αν, float τ_w, float refh, float h, float W)
        {
            Check_αν(αν);
            Check_τ(τ_w);
            Check_h(h);
            Check_h(refh);
            Check_W(W);


            return αν * (τ_w / (ρ * g)) * (W / refh) * Log(refh / h);

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="L_basin"></param>
        /// <param name="n"></param>
        /// <param name="h"></param>
        /// <returns></returns>
        public static float NaturalPeriodOscillationOpenBay(float L_basin, float n, float h)
        {
            Check_L(L_basin);
            Check_n(n);
            Check_n(h);

            return (4f * L_basin) / ((2f * n + 1f) * Sqrt(g * h));

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="L_basin"></param>
        /// <param name="n"></param>
        /// <param name="h"></param>
        /// <returns></returns>
        public static float NaturalPeriodOscillationClosedBasin(float L_basin, float n, float h)
        {
            Check_L(L_basin);
            Check_n(n);
            Check_n(h);

            return (2f * L_basin) / (n * Sqrt(g * h));

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="H"></param>
        /// <param name="δ"></param>
        /// <param name="K_D"></param>
        /// <param name="α"></param>
        /// <returns></returns>
        public static float HudsonsFormulaFor_D_n50(float H, float δ, float K_D, float α)
        {
            Check_H(H);
            Check_δ(δ);
            Check_K_D(K_D);
            Check_α(α);


            return H / Pow(δ * (K_D * α), (1f / 3f));

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="ρ_s"></param>
        /// <param name="D_n50"></param>
        /// <returns></returns>
        public static float WeightOfArmourUnit(float ρ_s, float D_n50)
        {
            Check_ρ_s(ρ_s);
            Check_D_n50(D_n50);

            return ρ_s * g * D_n50 * D_n50 * D_n50;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="ρ_s"></param>
        /// <returns></returns>
        public static float Find_DELTA(float ρ_s)
        {
            Check_ρ_s(ρ_s);

            return (ρ_s / ρ) - 1f;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="h"></param>
        /// <param name="γ_b"></param>
        /// <param name="H_BW"></param>
        /// <returns></returns>
        public static float BreakerCriteria(float h, float γ_b, float H_BW)
        {
            Check_h(h);
            Check_γ_b(γ_b);
            Check_H(H_BW);


            float H_b = h * γ_b;

            if (H_BW > H_b)
            {
                Debug.WriteLine("Wave is Broken");
            }
            else
            {
                Debug.WriteLine("Wave is NOT Broken");
                H_b = H_BW;
            }

            return H_b;


        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="Qstar"></param>
        /// <returns></returns>
        public static float OwensFormula_R_star(float a, float b, float Qstar)
        {

            Check_a(a);
            Check_b(b);

            return (-1f / b) * Log(Qstar / a);
            //R_star = (R_c/H_s) * Math.Sqrt(S_om/(2 * π)) * (1/GAMMA_y);

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="q"></param>
        /// <param name="H_s"></param>
        /// <param name="T_om"></param>
        /// <returns></returns>
        public static float OwensFormula_Q_star(float q, float H_s, float T_om)
        {
            Check_q(q);
            Check_H(H_s);
            Check_T_om(T_om);


            return q / (g * H_s * T_om);

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="H_stoe"></param>
        /// <param name="refL_om"></param>
        /// <returns></returns>
        public static float FindWaveSteepness(float H_stoe, float refL_om)
        {
            Check_H(H_stoe);
            Check_L(refL_om);


            return (H_stoe / refL_om);


        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="T_om"></param>
        /// <returns></returns>
        public static float Find_L_om(float T_om)
        {
            Check_T_om(T_om);

            return (g * T_om * T_om) / (2f * π);

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="R_star"></param>
        /// <param name="H_s"></param>
        /// <param name="γ_r"></param>
        /// <param name="S_om"></param>
        /// <returns></returns>
        public static float Find_R_c(float R_star, float H_s, float γ_r, float S_om)
        {
            Check_R_star(R_star);
            Check_H(H_s);
            Check_γ_r(γ_r);
            Check_S_om(S_om);

            return (R_star * H_s * γ_r) / Sqrt(S_om / (2f * π));

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="δ_z"></param>
        /// <param name="W"></param>
        /// <param name="B"></param>
        /// <param name="h_cl"></param>
        /// <returns></returns>
        public static float BruunsRule(float δ_z, float W, float B, float h_cl)
        {
            Check_δ_z(δ_z);
            Check_W(W);
            Check_B(B);
            Check_h_cl(h_cl);

            return δ_z * (W / (B + h_cl));
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="RCP_2_5_SLR_min"></param>
        /// <param name="RCP_2_5_SLR_max"></param>
        /// <returns></returns>
        public static float Find_δ_z_LowEmissions(float RCP_2_5_SLR_min, float RCP_2_5_SLR_max)
        {
            Check_RCP_2_5_SLR_min(RCP_2_5_SLR_min);
            Check_RCP_2_5_SLR_max(RCP_2_5_SLR_max);

            return (RCP_2_5_SLR_min + RCP_2_5_SLR_max) / 2;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="RCP_8_5_SLR_min"></param>
        /// <param name="RCP_8_5_SLR_max"></param>
        /// <returns></returns>
        public static float Find_DELTA_z_high_emissions(float RCP_8_5_SLR_min, float RCP_8_5_SLR_max)
        {
            Check_RCP_8_5_SLR_min(RCP_8_5_SLR_min);
            Check_RCP_8_5_SLR_max(RCP_8_5_SLR_max);

            return (RCP_8_5_SLR_min + RCP_8_5_SLR_max) / 2;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="K"></param>
        /// <param name="s"></param>
        /// <param name="H_b"></param>
        /// <param name="γ_b"></param>
        /// <param name="θ_b"></param>
        /// <returns></returns>
        public static float Find_Q_y(float K, float s, float H_b, float γ_b, float θ_b)
        {
            Check_K(K);
            Check_s(s);
            Check_H(H_b);
            Check_γ_b(γ_b);
            Check_θ_b(θ_b);


            return (K / (16f * (s - 1f) * Sqrt(γ_b))) * (Sqrt(g) * (Pow(H_b, (5f / 2f)))) * (Sind(2f * θ_b));

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="Q_grad"></param>
        /// <param name="p"></param>
        /// <param name="h_cl"></param>
        /// <param name="B"></param>
        /// <param name="SLR"></param>
        /// <returns></returns>
        public static float Find_Q_x(float Q_grad, float p, float h_cl, float B, float SLR)
        {
            Check_Q_grad(Q_grad);
            Check_p(p);
            Check_h_cl(h_cl);
            Check_B(B);
            Check_SLR(SLR);


            float grad_per_year = Q_grad * 365 * 24 * 60 * 60;
            return (grad_per_year) - (1 - p) * (h_cl + B) * (SLR);

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="Q_y"></param>
        /// <param name="refQ_y"></param>
        /// <param name="y"></param>
        /// <param name="refy"></param>
        /// <returns></returns>
        public static float Find_Q_grad(float Q_y, float refQ_y, float y, float refy)
        {
            Check_Q_y(Q_y);
            Check_Q_y(refQ_y);
            Check_y(y);
            Check_y(refy);

            return -(Q_y - refQ_y) / (y - refy);

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="Q_grad"></param>
        /// <param name="p"></param>
        /// <param name="h_cl"></param>
        /// <param name="B"></param>
        /// <param name="Q_x"></param>
        /// <returns></returns>
        public static float Find_SLR(float Q_grad, float p, float h_cl, float B, float Q_x)
        {
            Check_Q_grad(Q_grad);
            Check_p(p);
            Check_h_cl(h_cl);
            Check_B(B);
            Check_Q_x(Q_x);


            float grad_per_year = Q_grad * 365 * 24 * 60 * 60;
            return (Q_x + grad_per_year) / ((1 - p) * (h_cl + B));
        }

        public static float WaveLengthFromLambda(float ref_L, float T, float h)// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        {

            Check_L(ref_L);
            Check_T(T);
            Check_h(h);

            throw new NotImplementedException("Function is currently unavalible.");

            //return @(L_o)L_o - g * T ^ 2 / 2 / π * Math.Tanh(2 * π * h./ L_o);
            //L = fzero(public float_handle,ref_L);

        }

        #endregion

        #region Math Functions

        private static float Pow(float value, float power)
        {
            return (float)Math.Pow((double)value, (double)power);
        }

        private static float Exp(float value)
        {
            return (float)Math.Exp(value);
        }

        private static float Ln(float value)
        {
            return (float)Math.Log(value);
        }

        private static float Log(float value)
        {
            return (float)Math.Log10(value);
        }

        private static float Sqrt(float value)
        {
            return (float)Math.Sqrt(value);
        }

        private static float Cos(float value)
        {
            return (float)Math.Cos(value);
        }

        private static float Cosd(float value)
        {
            return (float)Math.Cos(value);//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        }

        private static float Cosh(float value)
        {
            return (float)Math.Cosh(value);
        }

        private static float Sin(float value)
        {
            return (float)Math.Sin(value);
        }

        private static float Sinh(float value)
        {
            return (float)Math.Sinh(value);
        }

        private static float Sind(float value) //<<<<<Convert <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        {
            return (float)Math.Sin(value);
        }

        private static float Asind(float value) //<<<<<Convert <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        {
            return (float)Math.Asin(value);
        }

        private static float Tan(float value)
        {
            return (float)Math.Tan(value);
        }

        private static float Tanh(float value)
        {
            return (float)Math.Tanh(value);
        }

        #endregion

        #region Check Values

        private static void Check_Q_x(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_Q_y(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_y(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_Q_grad(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_p(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_SLR(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_K(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_s(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_θ_b(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_RCP_8_5_SLR_max(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_RCP_8_5_SLR_min(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_RCP_2_5_SLR_max(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_RCP_2_5_SLR_min(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_δ_z(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_B(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_h_cl(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_R_star(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_γ_r(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_S_om(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_q(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_T_om(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_a(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_b(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_γ_b(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_ρ_s(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_D_n50(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_δ(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }
        private static void Check_K_D(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_α(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_αν(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_τ(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_W(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_C_10(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_U_10(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_ψ(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_R_n(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_z_100(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_γ(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_tanβ(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_β(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_K_r(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_K_s(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_θ(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_μ(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_ξ(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_σ(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_xtv(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_c(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
            if (value <= 0)
            {
                throw new ArithmeticException("Value is less than or equal to 0.");
            }
        }

        private static void Check_Ap(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_P_max(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_ps(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_K_p(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_η(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
        }

        private static void Check_N(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
            if (value <= 0)
            {
                throw new ArithmeticException("Value is less than or equal to 0.");
            }
        }

        private static void Check_n(float value) //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
            if (value <= 0)
            {
                throw new ArithmeticException("Value is less than or equal to 0.");
            }
        }

        private static void Check_H_rms(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
            if (value < 0)
            {
                throw new ArithmeticException("Value is less than 0.");
            }
        }

        private static void Check_H_avg(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
            if (value < 0)
            {
                throw new ArithmeticException("Value is less than 0.");
            }
        }

        private static void Check_H_s(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
            if (value < 0)
            {
                throw new ArithmeticException("Value is less than 0.");
            }
        }

        private static void Check_t(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
            if (value < 0)
            {
                throw new ArithmeticException("Value is less than 0.");
            }
        }

        private static void Check_x(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
            if (value < 0)
            {
                throw new ArithmeticException("Value is less than 0.");
            }
        }

        private static void Check_h(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
            if (value < 0)
            {
                throw new ArithmeticException("Value is less than 0.");
            }
        }

        private static void Check_z(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
            if (value > 0)
            {
                throw new ArithmeticException("Value is greater than 0.");
            }
        }

        private static void Check_k(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
            if (value <= 0)
            {
                throw new ArithmeticException("Value is less than or equal to 0.");
            }
        }

        private static void Check_H(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
            if (value <= 0)
            {
                throw new ArithmeticException("Value is less than or equal to 0.");
            }
        }

        private static void Check_ω(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
            if (value <= 0)
            {
                throw new ArithmeticException("Value is less than or equal to 0.");
            }
        }

        private static void Check_L(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
            if (value <= 0)
            {
                throw new ArithmeticException("Value is less than or equal to 0.");
            }
        }

        /// <summary>
        /// 
        /// [Symbol:  ]
        /// [Type  :  ]
        /// [Unit  :  ]
        /// </summary>
        /// <value>
        /// 
        /// </value>
        /// <remarks>
        /// 
        /// </remarks>    
        private static void Check_T(float value)
        {
            if (Object.ReferenceEquals(null, value))
            {
                throw new ArgumentNullException("Value is null.");
            }
            if (value <= 0)
            {
                throw new ArithmeticException("Value is less than or equal to 0.");
            }
        }

        #endregion

    }
}