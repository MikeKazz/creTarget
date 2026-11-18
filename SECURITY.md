# Security policy

## Supported versions

Security fixes are applied to the latest release and to `main`. Older releases
may receive a fix when the change is small and the release is still widely used;
they are not guaranteed support.

| Version | Supported |
|---|---|
| latest release | Yes |
| `main` | Yes |
| older releases | best effort |

## Reporting a vulnerability

Please do **not** open a public issue for a security problem.

Use GitHub's private vulnerability reporting for this repository:
[Security → Report a vulnerability](https://github.com/MikeKazz/creTarget/security/advisories/new).
If that form is unavailable, contact [@MikeKazz](https://github.com/MikeKazz)
privately through GitHub and ask for a secure channel.

Include, where possible:

- the affected version or commit,
- a minimal reproduction,
- the impact you believe is possible, and
- any suggested fix.

You can expect an acknowledgement within a few business days, an assessment,
and a coordinated disclosure once a fix is available. Please allow reasonable
time for a release before public disclosure.

## Scope

In scope:

- arbitrary code execution or command injection through pipeline inputs,
  parameters, file names, or downloaded assets;
- path traversal or unsafe archive extraction;
- bypass of the checksum or provenance checks for downloaded reference data,
  model checkpoints, or ATAC peaks;
- vulnerabilities in the container image or GitHub Actions workflows;
- credentials or tokens exposed by the codebase or CI configuration.

Out of scope:

- the biological performance of a guide, the accuracy of a published dataset,
  or the interpretation of a score;
- vulnerabilities in third-party tools (FlashFry, launch-dCas9, Nextflow,
  bedtools, samtools, Java) that are not caused by how creTarget invokes them;
- reports that require a compromised host, a malicious conda channel, or
  already-elevated privileges;
- missing hardening that has no demonstrated impact.

## Safe harbor

We will not pursue legal action against researchers who:

- make a good-faith effort to avoid privacy violations, data destruction, and
  service interruption;
- report the issue privately and give us a reasonable opportunity to fix it; and
- do not exploit the issue beyond what is needed to demonstrate it.
