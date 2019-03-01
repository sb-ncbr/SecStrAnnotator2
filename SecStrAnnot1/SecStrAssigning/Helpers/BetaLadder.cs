using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Cif.Components;

namespace protein.SecStrAssigning.Helpers
{
    class BetaLadder {
        public enum LadderType { Parallel, Antiparallel }
        public LadderType Type { get; private set; }

        public enum HBondDirection { From0To1, From1To0 }
        private HBondDirection InvertedDirection(HBondDirection d){
            return d == HBondDirection.From0To1 ? HBondDirection.From1To0 : HBondDirection.From0To1;
        }
        public HBondDirection FirstHBondDirection { get; set; }
        public HBondDirection LastHBondDirection { get; set; }

        public int Start0 { get; set; }
        public int End0 { get; set; }

        public int Start1 { get; set; }
        public int End1 { get; set; }

        // Z is defined as follows: Z=3*resi for backbone NH group, Z=3*resi+2 for backbone CO group
        public int ZStart0{ get { return 3 * Start0 + (FirstHBondDirection == HBondDirection.From1To0 ? 2 : 0); } }
        public int ZEnd0{ get { return 3 * End0 + (LastHBondDirection == HBondDirection.From1To0 ? 2 : 0); } }
        public int ZStart1 {
            get { 
                if (Type == LadderType.Parallel)
                    return 3 * Start1 + (FirstHBondDirection == HBondDirection.From0To1 ? 2 : 0);
                else
                    return 3 * Start1 + (LastHBondDirection == HBondDirection.From0To1 ? 2 : 0); 
            }
        }
        public int ZEnd1 {
            get { 
                if (Type == LadderType.Parallel)
                    return 3 * End1 + (LastHBondDirection == HBondDirection.From0To1 ? 2 : 0);
                else
                    return 3 * End1 + (FirstHBondDirection == HBondDirection.From0To1 ? 2 : 0); 
            }
        }

        private BetaLadder(LadderType type, int start0, int end0, int start1, int end1, HBondDirection firstHBondDirection, HBondDirection lastHBondDirection){
            Type=type;
            Start0 = start0;
            End0 = end0;
            Start1 = start1;
            End1 = end1;
            FirstHBondDirection = firstHBondDirection;
            LastHBondDirection = lastHBondDirection;
        }

        public BetaLadder(LadderType type, int res0, int res1, HBondDirection hBondDirection) 
            : this(type, res0, res0, res1, res1, hBondDirection, hBondDirection) {}

        public void AddOneHBond(){
            if (Type == LadderType.Antiparallel) {
                if (LastHBondDirection == HBondDirection.From0To1) {
                    LastHBondDirection = HBondDirection.From1To0;
                } else {
                    End0 += 2;
                    Start1 -= 2;
                    LastHBondDirection = HBondDirection.From0To1;
                }
            } else {
                if (LastHBondDirection == HBondDirection.From0To1) {
                    End1 += 2;
                    LastHBondDirection = HBondDirection.From1To0;
                } else {
                    End0 += 2;
                    LastHBondDirection = HBondDirection.From0To1;
                }
            }
        }

        public BetaLadder Inverted(){
            HBondDirection invFirstDir = (Type == LadderType.Antiparallel) ? InvertedDirection (LastHBondDirection) : InvertedDirection (FirstHBondDirection);
            HBondDirection invLastDir = (Type == LadderType.Antiparallel) ? InvertedDirection (FirstHBondDirection) : InvertedDirection (LastHBondDirection);
            return new BetaLadder (Type, Start1, End1, Start0, End0, invFirstDir,invLastDir);
        }

    }	
}

