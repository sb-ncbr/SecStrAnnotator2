using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using Cif;
using Cif.Tables;
using Cif.Filtering;
using Cif.Components;

namespace SecStrAnnotator2
{
    public class MyStopwatch {
        private DateTime t0;

        public MyStopwatch(){
            Start();
        }

        /** Reset to 0 */
        public void Start(){
            this.t0 = DateTime.Now;
        }

        /** Print the time since the last start/stop without reseting */
        public void Lap(string message){
            Lib.WriteLineDebug($"{message}: {DateTime.Now.Subtract(t0)}");
        }
        
        /** Print the time since the last start/stop and reset */
        public void Stop(string message){
            Lap(message);
            Start();
        }
    }
}