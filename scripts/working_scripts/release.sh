RELEASE='SecStrAnnotator_2.2'
FRAMEWORK='netcoreapp3.0'
PLUGIN_RELEASE='secstrapi_plugin-1.4-release.py'

DIR=`dirname $0`
cd $DIR/../..
echo $PWD

mkdir  $RELEASE
cp  bin/Release/$FRAMEWORK/SecStrAnnotator.dll  $RELEASE/
cp  bin/Release/netcoreapp3.0/SecStrAnnotator.runtimeconfig.json  $RELEASE/
cp  SecStrAnnotator_config-for_release.json  $RELEASE/SecStrAnnotator_config.json
cp  scripts/SecStrAnnotator_batch.py  $RELEASE/
cp  README.md  $RELEASE/
cp  LICENSE  $RELEASE/
cp  -r  examples  $RELEASE/

mkdir  $RELEASE/scripts
cp  scripts/SecStrAnnotator_batch.py  $RELEASE/scripts/
cp  scripts/SecStrAPI_pipeline.py  $RELEASE/scripts/
cp  scripts/SecStrAPI_pipeline_settings.json  $RELEASE/scripts/
cp  scripts/script_align.py  $RELEASE/scripts/
cp  scripts/script_session.py  $RELEASE/scripts/
cp  scripts/lib.py  $RELEASE/scripts/
cp  scripts/constants.py  $RELEASE/scripts/
cp  scripts/requirements.txt  $RELEASE/scripts/
cp  scripts/mkdssp  $RELEASE/scripts/
cp  -r  scripts/secstrapi_data_preparation  $RELEASE/scripts/
cp  -r  scripts/R_sec_str_anatomy_analysis  $RELEASE/scripts
rm  $RELEASE/scripts/R_sec_str_anatomy_analysis/*obsolete*
mkdir  $RELEASE/scripts/secstrapi_plugin
cp  scripts/secstrapi_plugin/$PLUGIN_RELEASE  $RELEASE/scripts/secstrapi_plugin

rm  -f  $RELEASE-release.zip
zip  $RELEASE-release.zip  -r $RELEASE
rm  -r  $RELEASE