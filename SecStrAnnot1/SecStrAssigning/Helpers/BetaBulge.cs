using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;

namespace protein.SecStrAssigning.Helpers
{
    class BetaBulge : IBulgeOrLadder {
        public enum BulgeType { Classic, Wide, Antiparallel22, Antiparallel33, Antiparallel15, Antiparallel23, Parallel14, Parallel32, Parallel13, Parallel33 } // New types might be added
        
        public BulgeType Type {get;private set;}
        public int StartShort{ get; private set; }
        public int EndShort{ get; private set; }
        public int StartLong{ get; private set; }
        public int EndLong{ get; private set; }
        public BetaBulge (BulgeType type, int startShort, int endShort, int startLong, int endLong){
            Type=type;
            StartShort=startShort;
            EndShort=endShort;
            StartLong=startLong;
            EndLong=endLong;
        }
    }	
}

