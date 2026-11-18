# Support

## Read this first

- [Documentation](https://mikekazz.github.io/creTarget/) — installation,
  inputs, outputs, parameters, cluster setup, and troubleshooting.
- [FAQ](docs/guides/faq.md) — short answers to common questions.
- [Troubleshooting](docs/guides/troubleshooting.md) — errors with a known cause.
- [CHANGELOG](CHANGELOG.md) — behavior changes between releases.

## Ask a question

Open a [question issue](https://github.com/MikeKazz/creTarget/issues/new?template=question.yml).
Include the version (`cretarget --version` and `nextflow -version`), the command
you ran, and what you expected to happen.

## Report a bug

Use the [bug report template](https://github.com/MikeKazz/creTarget/issues/new?template=bug_report.yml).
A useful report contains:

1. creTarget version, Nextflow version, Java version, and execution profile.
2. The command line, with local paths redacted if necessary.
3. `cretarget doctor` output.
4. The relevant part of `results/pipeline_info/trace.txt` and
   `results/run_manifest.json`.
5. The smallest input that reproduces the problem, or a description of how to
   construct one.

Please do not paste private data, credentials, or full genome-scale outputs into
an issue.

## Report a security issue

Follow [SECURITY.md](SECURITY.md). Do not use a public issue.

## Commercial or collaborative support

There is no commercial support offering. For collaborations, methods questions,
or benchmark requests, contact the maintainer through GitHub.

## Response expectations

This is a research tool maintained on a best-effort basis. Issues are triaged
when time allows; a clear reproduction and a small input move to the front of
the queue.
