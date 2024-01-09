using Documenter

struct CustomDeployConfig <: Documenter.DeployConfig
    repo::String
    event::String
    ref::String
    function CustomDeployConfig()
        github_repository = get(ENV, "GITHUB_REPOSITORY", "") # "JuliaQuantumControl/QuantumControlExamples.jl"
        github_event_name = get(ENV, "GITHUB_EVENT_NAME", "") # "push", "pull_request" or "cron" (?)
        github_ref = get(ENV, "GITHUB_REF", "") # "refs/heads/$(branchname)" for branch, "refs/tags/$(tagname)" for tags
        return new(github_repository, github_event_name, github_ref)
    end
end

function Documenter.post_status(::CustomDeployConfig; kwargs...)
    return Documenter.post_status(Documenter.GitHubActions(); kwargs...)
end

function Documenter.authenticated_repo_url(::CustomDeployConfig; kwargs...)
    return Documenter.authenticated_repo_url(Documenter.GitHubActions())
end

RX_DATE = r"20\d\d-\d\d$"

function Documenter.deploy_folder(cfg::CustomDeployConfig; kwargs...)
    all_ok = true
    if isempty(cfg.repo) || isempty(cfg.repo)
        @info "Do not deploy: Not running on Github Actions"
        all_ok = false
    end
    if !endswith(cfg.ref, RX_DATE)
        @info "Do not deploy: ref=$(repr(cfg.ref)) does not end with YYYY-MM"
        all_ok = false
    end
    all_ok && @info "OK to deploy"
    branch = "gh-pages"
    is_preview = false
    repo = "github.com/$(cfg.repo)"
    subfolder = split(cfg.ref, "/")[end]
    deploy_decision = Documenter.DeployDecision(all_ok, branch, is_preview, repo, subfolder)
    @show deploy_decision # DEBUG
    return deploy_decision
end

function get_git_versions(from="tag")
    cmd = `git $from --format="%(refname:short)"`
    output = read(cmd, String)
    return split(strip(output), '\n')
end


# Get all the versions that should be in the versions menu
function get_versions()
    folders = sort(
        filter(
            v -> startswith(v, RX_DATE),
            collect(Set([get_git_versions("tag")..., get_git_versions("branch")...]))
        ),
        rev=true
    )
    if isempty(folders)
        return nothing
    else
        return ["latest" => folders[begin], [(v => v) for v in folders]...]
    end
end



deploydocs(;
    repo="github.com/JuliaQuantumControl/QuantumControlExamples.jl",
    devbranch="master",
    deploy_config=CustomDeployConfig(),
    versions=get_versions(),
)
