"""prepareSlurmScriptExecution: generates an SBATCH array-job submit file.

Everything here is string generation plus two file reads, so the assertions
check the generated directives against Slurm's actual syntax rather than
against whatever the function currently emits.
"""
import pytest

from hera.utils.slurm import prepareSlurmScriptExecution


@pytest.fixture()
def jobs_file(tmp_path):
    """A job-directory list with three entries, so the array must be 1-3."""
    path = tmp_path / "jobPaths.txt"
    path.write_text("/jobs/a\n/jobs/b\n/jobs/c\n", encoding="utf-8")
    return path


@pytest.fixture()
def submit_path(tmp_path):
    return tmp_path / "submit_all.sh"


def _generate(jobs_file, submit_path, **kwargs):
    kwargs.setdefault("script", "echo hello")
    prepareSlurmScriptExecution(
        jobDirListFilePath=str(jobs_file),
        slurmExecutionFilePath=str(submit_path),
        **kwargs,
    )
    return submit_path.read_text(encoding="utf-8")


@pytest.mark.unit
class TestArgumentValidation:
    """Exactly one of script / scriptPath, and neither is silently accepted."""

    def test_neither_script_nor_path_writes_nothing(self, jobs_file, submit_path):
        prepareSlurmScriptExecution(
            jobDirListFilePath=str(jobs_file),
            slurmExecutionFilePath=str(submit_path),
        )
        assert not submit_path.exists()

    def test_both_script_and_path_writes_nothing(self, jobs_file, submit_path, tmp_path):
        script_file = tmp_path / "s.sh"
        script_file.write_text("echo hi", encoding="utf-8")
        prepareSlurmScriptExecution(
            script="echo hi",
            scriptPath=str(script_file),
            jobDirListFilePath=str(jobs_file),
            slurmExecutionFilePath=str(submit_path),
        )
        assert not submit_path.exists()

    def test_missing_job_list_writes_nothing(self, submit_path, tmp_path):
        prepareSlurmScriptExecution(
            script="echo hi",
            jobDirListFilePath=str(tmp_path / "absent.txt"),
            slurmExecutionFilePath=str(submit_path),
        )
        assert not submit_path.exists()

    def test_script_from_a_file_is_inlined(self, jobs_file, submit_path, tmp_path):
        script_file = tmp_path / "s.sh"
        script_file.write_text("echo from-file", encoding="utf-8")
        prepareSlurmScriptExecution(
            scriptPath=str(script_file),
            jobDirListFilePath=str(jobs_file),
            slurmExecutionFilePath=str(submit_path),
        )
        assert "echo from-file" in submit_path.read_text(encoding="utf-8")


@pytest.mark.unit
class TestGeneratedDirectives:
    def test_shebang_comes_first(self, jobs_file, submit_path):
        assert _generate(jobs_file, submit_path).startswith("#!/bin/bash")

    def test_array_size_matches_the_job_count(self, jobs_file, submit_path):
        """Three job directories must produce a 1-3 array, not 1-2 or 1-4."""
        assert "#SBATCH --array=1-3" in _generate(jobs_file, submit_path)

    def test_array_size_follows_the_file(self, tmp_path, submit_path):
        jobs = tmp_path / "one.txt"
        jobs.write_text("/only/one\n", encoding="utf-8")
        assert "#SBATCH --array=1-1" in _generate(jobs, submit_path)

    def test_job_name_is_used(self, jobs_file, submit_path):
        text = _generate(jobs_file, submit_path, jobName="myRun")
        assert "#SBATCH --job-name=myRun" in text

    def test_default_job_name(self, jobs_file, submit_path):
        assert "#SBATCH --job-name=generic_job" in _generate(jobs_file, submit_path)

    def test_the_script_body_is_appended(self, jobs_file, submit_path):
        assert "echo hello" in _generate(jobs_file, submit_path)

    def test_the_task_id_selects_a_directory(self, jobs_file, submit_path):
        """Each array task must cd into its own directory, by line number."""
        text = _generate(jobs_file, submit_path)
        assert 'sed -n "${SLURM_ARRAY_TASK_ID}p"' in text
        assert str(jobs_file) in text


@pytest.mark.unit
class TestOptionalDirectives:
    def test_processors_directive_appears_only_when_asked(self, jobs_file, submit_path):
        assert "#SBATCH -n 4" in _generate(
            jobs_file, submit_path, allocateProcessorsPerRun=4
        )
        assert "#SBATCH -n" not in _generate(jobs_file, submit_path)

    def test_exclusive_directive_appears_only_when_asked(self, jobs_file, submit_path):
        assert "#SBATCH --exclusive" in _generate(jobs_file, submit_path, exclusive=True)
        assert "--exclusive" not in _generate(jobs_file, submit_path, exclusive=False)

    def test_quiet_suppresses_the_progress_echo(self, jobs_file, submit_path):
        assert "Running script for job" in _generate(jobs_file, submit_path)
        assert "Running script for job" not in _generate(
            jobs_file, submit_path, quiet=True
        )

    def test_memory_directive_appears_only_when_asked(self, jobs_file, submit_path):
        assert "mem=" not in _generate(jobs_file, submit_path)

    @pytest.mark.xfail(
        strict=True,
        reason="B7: the memory directive is emitted as '#SBATCH -mem=8G'; Slurm's "
               "option is '--mem'. A single dash makes sbatch reject the script. "
               "See the consolidated findings issue.",
    )
    def test_memory_directive_uses_slurms_option_name(self, jobs_file, submit_path):
        text = _generate(jobs_file, submit_path, memoryInGB=8)
        assert "#SBATCH --mem=8G" in text
