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
    class MOMAnnotationContext
    {
        public AnnotationContext Context { get; private set; }
        public (int, int, int)[] TemplateFST { get; private set; } // tuple<First,Second,Type>
        public (int, int, int)[] CandidateFST { get; private set; }
        public MOMAnnotationContext(AnnotationContext inner)
        {
            inner.ValidateOrdering();
            Context = inner;
            List<(int, int, int)> templateConnections = new List<(int, int, int)>();
            for (int i = 0; i < Context.Templates.Length; i++)
            {
                for (int j = i + 1; j < Context.Templates.Length; j++)
                {
                    if (Context.TemplateConnectivity[i, j] != 0)
                    {
                        templateConnections.Add((i, j, Context.TemplateConnectivity[i, j]));
                    }
                }
            }
            List<(int, int, int)> candidateConnections = new List<(int, int, int)>();
            for (int i = 0; i < Context.Candidates.Length; i++)
            {
                for (int j = i + 1; j < Context.Candidates.Length; j++)
                {
                    if (Context.CandidateConnectivity[i, j] != 0)
                    {
                        candidateConnections.Add((i, j, Context.CandidateConnectivity[i, j]));
                    }
                }
            }

            TemplateFST = Context.Templates.IndicesWhere(sse => sse.IsHelix).Select(i => (i, i, 0)) //helix nodes
                .Concat(templateConnections.Select(t => (Math.Min(t.Item1, t.Item2), Math.Max(t.Item1, t.Item2), t.Item3))) //ladder nodes
                .ToArray();
            CandidateFST = Context.Candidates.IndicesWhere(sse => sse.IsHelix).Select(i => (i, i, 0)) //helix nodes
                .Concat(candidateConnections.Select(t => (Math.Min(t.Item1, t.Item2), Math.Max(t.Item1, t.Item2), t.Item3))) //ladder nodes
                .ToArray();
        }
    }
}