using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

using Cif.Components;
using protein.Sses;
using protein.Libraries;

namespace protein.Annotating {
    /** Run primary annotator; if it runs longer than timeoutSeconds, show warning and run backup annotator. */
    class FallbackAnnotator : IAnnotator {
        public ICancellableAnnotator PrimaryAnnotator;
        public double TimeoutSeconds;
        public string Warning;
        public IAnnotator BackupAnnotator;

        public FallbackAnnotator(ICancellableAnnotator primaryAnnotator, double timeoutSeconds, string timeoutWarning, IAnnotator backupAnnotator) {
            PrimaryAnnotator = primaryAnnotator;
            TimeoutSeconds = timeoutSeconds;
            Warning = timeoutWarning;
            BackupAnnotator = backupAnnotator;
        }

        public List<(int, int)> GetMatching() {
            using(var tokenSource = new CancellationTokenSource()){
                CancellationToken cancellationToken = tokenSource.Token;
                var task = Task.Run(() => PrimaryAnnotator.GetMatching(cancellationToken), cancellationToken);
                tokenSource.CancelAfter(TimeSpan.FromSeconds(TimeoutSeconds));
                try {
                    task.Wait();
                    return task.Result;
                } catch (Exception ex) when (ex is OperationCanceledException || ex is AggregateException) {
                    if (task.Status == TaskStatus.Canceled) {
                        // OK
                        if (Warning != null) {
                            Lib.WriteWarning(Warning);
                        }
                        return BackupAnnotator.GetMatching();
                    } else {
                        throw;
                    }
                }
            }
        }

    }
}