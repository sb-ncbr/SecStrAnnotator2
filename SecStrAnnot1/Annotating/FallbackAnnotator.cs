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
    /** Run primary annotator; if it runs longer than timeoutSeconds, show warning and run backup annotator. */
    class FallbackAnnotator : IAnnotator
    {
        public IAnnotator PrimaryAnnotator;
        public double TimeoutSeconds;
        public string Warning;
        public IAnnotator BackupAnnotator;

        public FallbackAnnotator(IAnnotator primaryAnnotator, double timeoutSeconds, string timeoutWarning, IAnnotator backupAnnotator)
        {
            PrimaryAnnotator = primaryAnnotator;
            TimeoutSeconds = timeoutSeconds;
            Warning = timeoutWarning;
            BackupAnnotator = backupAnnotator;
        }

        public List<(int, int)> GetMatching()
        {
            var task = Task.Run(() => PrimaryAnnotator.GetMatching());
            if (task.Wait(TimeSpan.FromSeconds(TimeoutSeconds)))
            {
                return task.Result;
            }
            else
            {
                if (Warning != null)
                {
                    Lib.WriteWarning(Warning);
                }
                return BackupAnnotator.GetMatching();
            }
        }
    }
}