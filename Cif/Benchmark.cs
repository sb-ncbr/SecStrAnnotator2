using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using Cif.Libraries;

namespace /*SecStrAnnot2.*/Cif
{
    public class Benchmark
    {
        public const int N = 40000;

        public class One : Interface {
            private int[] array;
            public One() { array = Enumerable.Range(0, N).ToArray(); }
            // [MethodImpl(MethodImplOptions.NoInlining)]
            public int Get(int i){ return array[i]; }
        }

        public class Super : Interface {
            private int[] array;
            public Super() { array = Enumerable.Range(0, N).ToArray(); }
            // [MethodImpl(MethodImplOptions.NoInlining)]
            public virtual int Get(int i){ return array[i]; }
        }

        public class Sub : Super {
            private int[] array;
            public Sub() { array = Enumerable.Range(0, N).ToArray(); }
            // [MethodImpl(MethodImplOptions.NoInlining)]
            // public override int Get(int i){ return array[i]; }
        }

        public interface Interface {
            int Get(int i);
        }

        public static void RunOn(Interface instance){
            instance.Get(0);
            DateTime t0 = DateTime.Now;
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    int x = instance.Get(j);
                }
            }
            DateTime t1 = DateTime.Now;
            Func<TimeSpan,string> Format = span => span.TotalSeconds.ToString("0.000");
            Lib.WriteLineDebug(" " + Format(t1-t0) + "  " + instance.GetType());
        }

        public static void RunOnOne(One instance){
            instance.Get(0);
            DateTime t0 = DateTime.Now;
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    int x = instance.Get(j);
                }
            }
            DateTime t1 = DateTime.Now;
            Func<TimeSpan,string> Format = span => span.TotalSeconds.ToString("0.000");
            Lib.WriteLineDebug(" " + Format(t1-t0) + "  " + instance.GetType());
        }

        public static void RunOnSuper(Super instance){
            instance.Get(0);
            DateTime t0 = DateTime.Now;
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    int x = instance.Get(j);
                }
            }
            DateTime t1 = DateTime.Now;
            Func<TimeSpan,string> Format = span => span.TotalSeconds.ToString("0.000");
            Lib.WriteLineDebug(" " + Format(t1-t0) + "  " + instance.GetType());
        }

        public static void Run() {
            One one = new One();
            Super super = new Super();
            Sub sub = new Sub();
            Super incognito = new Sub();

            RunOnOne(one);
            RunOn(one);
            RunOn(super);
            RunOn(sub);
            RunOn(incognito);
            RunOnSuper(super);
            RunOnSuper(sub);
        }
    }
}