# Releasing `genome-sampler`

1. Create a new git branch, `r20WX.YZ`
2. Create a new commit on the new branch (`-m 'REL: 20WX.YZ'`), with the
   following changes:
    1. Update `.github/workflows/main.yml`:
        1. Replacing all instances of `master` with the release branch name
           (`r20WX.YZ`)
        2. Set `build-target: release` in `build-and-test` job
        3. Update `envFile` to point to desired release in `integrate` job
        4. Update env file URL in `integrate` job
        5. Update `-c qiime2` in `integrate` job
        6. Update `if github.ref == 'refs/heads/r20WX.YZ` in `doc-publish` job
    2. Update `docs/tutorial.md`:
        1. Replace all instances of the prior release version string with
           `20WX.YZ`
3. Tag the commit created in step 2, `git tag 20WX.YZ`
4. Push the new commit and tag up to
   https://github.com/caporaso-lab/genome-sampler

## diff format of the changes above

    diff --git a/.github/workflows/main.yml b/.github/workflows/main.yml
    index 35b3b32..62790af 100644
    --- a/.github/workflows/main.yml
    +++ b/.github/workflows/main.yml
    @@ -3,12 +3,12 @@ name: ci
     on:
       push:
         branches:
    -      - master
    +      - r2020.8
         paths-ignore:
           - 'conda-env-files/**'
       pull_request:
         branches:
    -      - master
    +      - r2020.8

     jobs:
       lint:
    @@ -53,7 +53,7 @@ jobs:

       integrate:
         needs: build-and-test
    -    if: github.ref == 'refs/heads/master'
    +    if: github.ref == 'refs/heads/r2020.8'
         strategy:
           matrix:
             os: [ubuntu-latest, macos-latest]
    @@ -79,7 +79,7 @@ jobs:
             run: |
               # TODO: once we catch up to the QIIME 2 release schedule, fix this
               envFile=qiime2-latest-py36-$QIIME2_PLATFORM-conda.yml
    -          wget https://raw.githubusercontent.com/qiime2/environment-files/master/latest/staging/$envFile
    +          wget https://data.qiime2.org/distro/core/$envFile
               sudo conda env create -q -p ./genome-sampler-env --file $envFile

           - name: install genome-sampler
    @@ -87,7 +87,7 @@ jobs:
               source "$CONDA/etc/profile.d/conda.sh"
               sudo conda install -p ./genome-sampler-env -q \
                 -c https://packages.qiime2.org/qiime2/latest/tested/ \
    +            -c qiime2 \
                 -c conda-forge \
                 -c bioconda \
                 -c defaults \
    @@ -111,7 +111,7 @@ jobs:

       commit-env-files:
         needs: integrate
    -    if: github.ref == 'refs/heads/master'
    +    if: github.ref == 'refs/heads/r2020.8'
         runs-on: ubuntu-latest
         steps:
           - uses: actions/checkout@v2
    @@ -181,7 +181,7 @@ jobs:
       doc-publish:
         needs: doc-build
         runs-on: ubuntu-latest
    -    if: github.ref == 'refs/heads/some-release-branch'
    +    if: github.ref == 'refs/heads/r2020.8'
         steps:
         - name: fetch built docs
           uses: actions/download-artifact@v1
    diff --git a/docs/tutorial.md b/docs/tutorial.md
    index 7aea345..594262a 100644
    --- a/docs/tutorial.md
    +++ b/docs/tutorial.md
    @@ -25,7 +25,7 @@ For linux installation environments, please run:

     ```bash
     wget https://raw.githubusercontent.com/caporaso-lab/genome-sampler/some-release-branch/conda-env-files/genome-sampler-py36-linux-conda.yml
    -conda env create -n genome-sampler-2020.6 --file genome-sampler-py36-linux-conda.yml
    +conda env create -n genome-sampler-2020.8 --file genome-sampler-py36-linux-conda.yml
     rm genome-sampler-py36-linux-conda.yml
     ```

    @@ -33,14 +33,14 @@ For macOS installation environments, please run:

     ```bash
     wget https://raw.githubusercontent.com/caporaso-lab/genome-sampler/some-release-branch/conda-env-files/genome-sampler-py36-osx-conda.yml
    -conda env create -n genome-sampler-2020.6 --file genome-sampler-py36-osx-conda.yml
    +conda env create -n genome-sampler-2020.8 --file genome-sampler-py36-osx-conda.yml
     rm genome-sampler-py36-osx-conda.yml
     ```

     ### Activate the `genome-sampler` conda environment

     ```bash
    -conda activate genome-sampler-2020.6
    +conda activate genome-sampler-2020.8
     ```

     ## Usage instructions

# Preparing a new development cycle

This is similar to what we do in "official" QIIME 2 packages - create an empty
commit with the version info, and tag it. This should be done one `master`.

```bash
git pull caporaso-lab master --tags
git commit --allow-empty -m 'VER: 20XY.Z.0.dev0'
git tag 20XY.Z.0.dev0
git push caporaso-lab master --tags
```
