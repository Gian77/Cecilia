// Nextflow does not coerce `--flag false` / `--flag true` CLI overrides to a
// real Boolean — they arrive as the Strings "false"/"true", and Groovy treats
// any non-empty String (including "false") as truthy. Every boolean param must
// be passed through asBool() before an `if (params.x)`-style check, or
// `--x false` silently takes the true branch.
class Utils {
    static boolean asBool(v) {
        return (v instanceof String) ? v.toBoolean() : (v as boolean)
    }

    // A process output declared `path("*.fastq")` emits a bare Path when
    // exactly one file matches, and a List<Path> when more than one does.
    // Indexing a bare Path with [n] silently returns its n-th path *component*
    // (Path implements Iterable<Path>) instead of erroring — wrap any [n]
    // access on such an output in asList() first to get a real element.
    static List asList(v) {
        return (v instanceof List) ? v : [v]
    }

    // -------------------------------------------------------------------------
    // USEARCH v11.0.667 was donated to the public domain under CC0-1.0
    // (rcedgar/usearch_old_binaries), so the pipeline fetches it on first run
    // instead of requiring every user to obtain a copy and set a host path.
    //
    // Why this runs from the entry workflow and not from a Nextflow process:
    // every USEARCH module interpolates `${params.usearch}` directly into its
    // script, and params cannot be rewritten at runtime — assigning
    // `params.usearch` from a workflow body is silently ignored (a no-op, not an
    // error; verified on Nextflow 26.04.3). So the binary must already exist at
    // the configured path before the first task is submitted. Threading it
    // through a channel instead would mean rewriting the input block of all
    // fifteen USEARCH modules and invalidating the whole -resume cache.
    // -------------------------------------------------------------------------
    static String ensureUsearch(bin, url, expectedSha, boolean fetch) {
        def target = new File(bin.toString())
        def sha    = expectedSha?.toString()

        // A binary already in place is never overwritten: it may be a build the
        // user pointed at deliberately via --usearch. Report a hash mismatch and
        // continue rather than second-guessing them.
        if (target.exists() && target.length() > 0) {
            if (sha) {
                def actual = sha256(target)
                if (actual != sha) {
                    System.err.println(
                        "WARNING: ${target} is not the USEARCH build this pipeline pins.\n" +
                        "         expected sha256 ${sha}\n" +
                        "         found    sha256 ${actual}\n" +
                        "         Continuing with it anyway.")
                }
            }
            if (!target.canExecute()) target.setExecutable(true, false)
            return target.absolutePath
        }

        def manual = "Fetch it by hand on a node with internet access, then re-run:\n" +
                     "    mkdir -p ${target.parent}\n" +
                     "    curl -fLo ${target} ${url}\n" +
                     "    chmod +x ${target}"

        if (!fetch) {
            throw new IllegalStateException(
                "USEARCH not found at ${target} and --usearch_fetch is false.\n${manual}")
        }

        target.parentFile?.mkdirs()

        // Download to a scratch name and only move it into place after the hash
        // checks out, so an interrupted or corrupted fetch can never leave a
        // truncated binary at the configured path for the next run to trust.
        def part = new File(target.parentFile,
                            "${target.name}.part.${UUID.randomUUID().toString().take(8)}")
        System.err.println("Fetching USEARCH v11 (CC0-1.0) -> ${target}")
        try {
            def conn = new URL(url.toString()).openConnection() as HttpURLConnection
            conn.setInstanceFollowRedirects(true)
            conn.setConnectTimeout(30000)
            conn.setReadTimeout(300000)
            conn.setRequestProperty('User-Agent', 'cecilia-nf')
            part.withOutputStream { out -> conn.inputStream.withStream { src -> out << src } }
        }
        catch (Exception e) {
            part.delete()
            throw new IllegalStateException(
                "Could not download USEARCH from ${url}\n  ${e.message}\n${manual}")
        }

        if (sha) {
            def got = sha256(part)
            if (got != sha) {
                part.delete()
                throw new IllegalStateException(
                    "Downloaded USEARCH failed its checksum — refusing to use it.\n" +
                    "  expected sha256 ${sha}\n  got      sha256 ${got}\n  from ${url}")
            }
        }

        part.setExecutable(true, false)
        java.nio.file.Files.move(part.toPath(), target.toPath(),
                                 java.nio.file.StandardCopyOption.REPLACE_EXISTING)
        System.err.println("USEARCH ready: ${target} (${target.length()} bytes, sha256 verified)")
        return target.absolutePath
    }

    private static String sha256(File f) {
        def md = java.security.MessageDigest.getInstance('SHA-256')
        f.withInputStream { s ->
            byte[] buf = new byte[65536]
            int n
            while ((n = s.read(buf)) > 0) md.update(buf, 0, n)
        }
        return md.digest().encodeHex().toString()
    }
}
