# GitHub Copilot instructions for PGL

Use original CSV records for atlas questions. Read [the shared access guide](../AGENTS.md) before retrieving catalogs or joining sampling information to a PGL result.

Prefer direct reads from the connected repository or checkout. For web browsing, preserve `?plain=1` on GitHub CSV URLs; the guide lists explicit alternative CSV endpoints. Verify that records were received, not merely a header or displayed line count. If every available method fails, report the exact failed URLs and request the required CSV or connected checkout.

Join result `structurename` to catalog `structure_name` exactly, preserving hemisphere and subdivisions. `sample_count` counts sampled wells; `donor_count` counts donor brains with samples; result `brainFreq` counts donors where PGL selected the structure for that gene list. Keep those measures separate. Cite the source CSV records and identify unmatched or ambiguous keys explicitly.

These instructions are for GitHub Copilot features that support repository instructions. Do not assume a separate Microsoft Copilot web chat loads this file when a GitHub URL is mentioned.
