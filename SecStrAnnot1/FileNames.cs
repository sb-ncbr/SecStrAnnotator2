using System.IO;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

using protein.Geometry;

namespace protein
{
    public class FileNames
    {
        public readonly string TemplateStructure;
        public readonly string TemplateAnnotatedHelices;

        public readonly string QueryStructure;
        public readonly string QueryAlignedStructure;
        public readonly string QueryAlignment;
        public readonly string QueryRenumberedStructure;
        public readonly string QueryDSSP;
        public readonly string QueryInputHelices;
        public readonly string QueryDetectedHelices;
        public readonly string QueryJoinedHelices;
        public readonly string QueryAnnotatedHelices;
        public readonly string QueryRmsds;
        public readonly string QueryLabel2Auth;
        public readonly string QuerySequence;

        public FileNames(string directory, string templateID, string queryID)
        {
            if (templateID != null){
                TemplateStructure = Path.Combine(Setting.Directory, templateID + Setting.PDB_FILE_EXT);
                TemplateAnnotatedHelices = Path.Combine(Setting.Directory, templateID + Setting.TEMPLATE_ANNOTATION_FILE_EXT);
            }

            if (queryID != null){
                QueryStructure = Path.Combine(Setting.Directory, queryID + Setting.PDB_FILE_EXT);
                QueryAlignedStructure = Path.Combine(Setting.Directory, queryID + Setting.ALIGNED_PDB_FILE_EXT);
                QueryAlignment = Path.Combine(Setting.Directory, queryID + Setting.ALIGNMENT_FILE_EXT);
                QueryRenumberedStructure = Path.Combine(Setting.Directory, queryID + Setting.RENUMBERED_PDB_FILE_EXT);
                QueryDSSP = Path.Combine(Setting.Directory, queryID + Setting.DSSP_OUTPUT_FILE_EXT);
                QueryInputHelices = Path.Combine(Setting.Directory, queryID + Setting.INPUT_SSES_FILE_EXT);
                QueryDetectedHelices = Path.Combine(Setting.Directory, queryID + Setting.DETECTED_SSES_FILE_EXT);
                QueryJoinedHelices = Path.Combine(Setting.Directory, queryID + Setting.JOINED_SSES_FILE_EXT);
                QueryAnnotatedHelices = Path.Combine(Setting.Directory, queryID + Setting.ANNOTATION_FILE_EXT);
                QueryRmsds = Path.Combine(Setting.Directory, queryID + Setting.RMSDS_FILE_EXT);
                QueryLabel2Auth = Path.Combine(Setting.Directory, queryID + Setting.LABEL2AUTH_FILE_EXT);
                QuerySequence = Path.Combine(Setting.Directory, queryID + Setting.SEQUENCE_FILE_EXT);
            }
        }

    }
}

