# Release and archival procedure

This procedure keeps a software tag, container digest, and Zenodo record tied to the same clean Git revision.

## Pre-release checks

1. Confirm that `uniscflow.__version__`, `CITATION.cff`, `.zenodo.json`, `CHANGELOG.md`, and the proposed `vX.Y.Z` tag agree.
2. Run the full unit suite and shell syntax checks.
3. Build from a clean checkout with `bash docker/build_current_image.sh`.
4. Verify the image independently with `bash docker/verify_image.sh IMAGE EXPECTED_COMMIT EXPECTED_VERSION`.
5. Confirm that the release repository and its GitHub Actions have the required container-package permissions. Use one Zenodo archival method for each release; do not create both a manual deposit and an automatic GitHub deposit for the same release.

## Release

Choose a new, unused version and update the source and metadata consistently. Create and push an annotated tag only after the checks pass:

```bash
VERSION=X.Y.Z  # Replace with the new version recorded in the source and metadata.
git tag -a "v${VERSION}" -m "UniScFlow ${VERSION}"
git push origin "v${VERSION}"
```

The release workflow verifies the source, builds the Python distributions, publishes the Linux/AMD64 image to GitHub Container Registry with SBOM and provenance attestations, and creates the GitHub release. It does not publish a Zenodo deposit.

Version 1.0.0 is archived at `10.5281/zenodo.22271248`; that DOI does not identify later changes on `main`. For a new release, create a new Zenodo version and upload its source archive, wheel, source provenance, and checksums from the clean release commit. Confirm the commit and repository URL, and download the uploaded files again to verify checksums before publishing. Update citation metadata for the new archive, and do not also enable automatic GitHub archiving for that same release.

## Post-release checks

1. Record the Git tag object, commit, GitHub release URL, GHCR digest, and Zenodo version DOI.
2. Confirm that the repository, release assets, container, and Zenodo record are accessible without authentication.
3. Verify that the version DOI agrees with the repository citation metadata. Make subsequent changes on `main`; do not move or recreate a released tag.

`latest` is a convenience tag and is not an archival identifier. The version tag, Git commit, container digest, and Zenodo DOI are the frozen identifiers.
