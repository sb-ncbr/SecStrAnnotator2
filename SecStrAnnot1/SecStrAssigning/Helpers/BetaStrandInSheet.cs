using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;

namespace protein.SecStrAssigning.Helpers
{
    class BetaStrandInSheet {
        public SSE SSE { get; set; }
        public int SheetId { get; set; }
        //public int Level { get; set; }
        public bool EvenUp { get; set; }
        public List<BetaStrandInSheet> UpNeighbours { get; set; }
        public List<BetaStrandInSheet> DownNeighbours { get; set; }
        // The ladders must be stored so that strand 0 belongs to this!
        //TODO a check of this condition should be included in AddUpLadder/AddDownLadder
        public List<BetaLadder> UpLadders { get; set; } 
        public List<BetaLadder> DownLadders{ get; set; }
        private bool DFSFlag { get; set; }

        public BetaStrandInSheet(SSE sse, int sheetId, /*int level, */bool evenUp){
            SSE=sse;
            SheetId=sheetId;
            //Level=level;
            EvenUp=evenUp;
            UpNeighbours=new List<BetaStrandInSheet>();
            DownNeighbours=new List<BetaStrandInSheet>();
            UpLadders=new List<BetaLadder>();
            DownLadders=new List<BetaLadder>();
            DFSFlag=false;
        }

        private void DFSAux(bool flagDiscovered, Action<BetaStrandInSheet> actionOnDiscover, Action<BetaStrandInSheet> actionOnReturn){
            if (DFSFlag != flagDiscovered) {
                DFSFlag = flagDiscovered;
                //if (flagDiscovered) Lib.WriteDebug ("[{0}], ", this.SSE.Label);
                actionOnDiscover (this);
                foreach (BetaStrandInSheet s in UpNeighbours.Union (DownNeighbours))
                    s.DFSAux (flagDiscovered, actionOnDiscover, actionOnReturn);
                actionOnReturn (this);
            }
        }

        public void DFS (Action<BetaStrandInSheet> actionOnDiscover, Action<BetaStrandInSheet> actionOnReturn){
                    //Lib.WriteLineDebug ("DFS from {0}:", this.SSE.Label);
            DFSAux (true, actionOnDiscover, actionOnReturn);
            DFSAux (false, s => {}, s => {});
            //Lib.WriteLineDebug ("");
        }

        public void DFS (Action<BetaStrandInSheet> actionOnDiscover){
            DFS (actionOnDiscover, s => {});
        }

        public override String ToString(){
            return SSE.Label + "[" + SSE.ChainID + SSE.Start + "-" + SSE.End + "]";
        }
    }
}

