using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Threading.Tasks;

using Cif.Components;
using protein.Sses;
using protein.Libraries;

namespace protein.Annotating
{
    class NiceAnnotatorWrapperWithCorrections : NiceAnnotatorWrapper
    {
        private List<List<String>> corrs; // each list contains [pdbid, label, chainid, startresi, endresi]
        Protein protein;

        public NiceAnnotatorWrapperWithCorrections(
            AnnotationContext context,
            Func<AnnotationContext, IAnnotator> innerAnnotatorConstructor,
            String correctionsFile,
            String pdbId,
            Protein p)
            : base(context, innerAnnotatorConstructor)
        {
            using (StreamReader r = new StreamReader(correctionsFile))
            {
                this.corrs = Lib.ReadCSV(r, 5, '\t', '#').Where(l => l[0] == pdbId).ToList();
            }
            this.protein = p;
        }
        
        public override IEnumerable<SseInSpace> GetAnnotatedCandidates()
        {
            SseInSpace[] result = base.GetAnnotatedCandidates().ToArray();
            for (int i = 0; i < result.Length; i++)
            {
                List<String> correction = corrs.Where(l => l[1] == Context.Templates[i].Label).DefaultIfEmpty(null).FirstOrDefault();
                if (correction != null)
                {
                    SseInSpace template = Context.Templates[i];
                    string chainID = correction[2];
                    int start = Int32.Parse(correction[3]);
                    int end = Int32.Parse(correction[4]);
                    Sse correctedSSE = new Sse(template.Label, chainID, start, end, template.Type, null);
                    List<double>[] dump;
                    result[i] = (start == 0 && end == 0) ?
                        SseInSpace.NewNotFound(template.Label)
                        : LibAnnotation.SSEsAsLineSegments_GeomVersion(protein.GetChain(correctedSSE.ChainID), new Sse[] { correctedSSE }.ToList(), out dump).First();
                }
            }
            return result;
        }
    }
}