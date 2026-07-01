//! Optional GUI frame profiling hooks.
//!
//! Normal GENtle builds keep these helpers inert. Builds with
//! `--features gui-profiler` can set `GENTLE_GUI_PROFILE=1` to stream Puffin
//! profile data over localhost for inspection with `puffin_viewer`.

use std::env;

pub const GUI_PROFILER_ENV: &str = "GENTLE_GUI_PROFILE";
pub const GUI_PROFILER_ADDR_ENV: &str = "GENTLE_GUI_PROFILE_ADDR";
pub const DEFAULT_GUI_PROFILER_ADDR: &str = "127.0.0.1:8585";

#[cfg(feature = "gui-profiler")]
static GUI_PROFILER_SERVER: std::sync::OnceLock<Option<puffin_http::Server>> =
    std::sync::OnceLock::new();

pub(crate) fn env_flag_value_enabled(value: Option<&str>) -> bool {
    value
        .map(|value| {
            matches!(
                value.trim().to_ascii_lowercase().as_str(),
                "1" | "true" | "yes" | "on"
            )
        })
        .unwrap_or(false)
}

pub fn init_from_env() {
    if !env_flag_value_enabled(env::var(GUI_PROFILER_ENV).ok().as_deref()) {
        return;
    }
    init_enabled_from_env();
}

#[cfg(feature = "gui-profiler")]
fn init_enabled_from_env() {
    let addr = env::var(GUI_PROFILER_ADDR_ENV)
        .ok()
        .filter(|addr| !addr.trim().is_empty())
        .unwrap_or_else(|| DEFAULT_GUI_PROFILER_ADDR.to_string());
    let _ = GUI_PROFILER_SERVER.get_or_init(|| {
        puffin::set_scopes_on(true);
        match puffin_http::Server::new(&addr) {
            Ok(server) => {
                eprintln!(
                    "I GENtle GUI profiler: streaming Puffin data on {addr}; run puffin_viewer to connect"
                );
                Some(server)
            }
            Err(err) => {
                eprintln!(
                    "W GENtle GUI profiler: could not listen on {addr}: {err}; profiling scopes remain enabled without HTTP streaming"
                );
                None
            }
        }
    });
}

#[cfg(not(feature = "gui-profiler"))]
fn init_enabled_from_env() {
    eprintln!(
        "W GENtle GUI profiler: {GUI_PROFILER_ENV}=1 ignored; rebuild with --features gui-profiler"
    );
}

pub fn begin_frame() {
    #[cfg(feature = "gui-profiler")]
    {
        if puffin::are_scopes_on() {
            puffin::GlobalProfiler::lock().new_frame();
        }
    }
}

#[macro_export]
macro_rules! gentle_gui_profile_scope {
    ($name:expr) => {
        #[cfg(feature = "gui-profiler")]
        ::puffin::profile_scope!($name);
    };
    ($name:expr, $data:expr) => {
        #[cfg(feature = "gui-profiler")]
        ::puffin::profile_scope!($name, $data);
    };
}

#[cfg(test)]
mod tests {
    use super::env_flag_value_enabled;

    #[test]
    fn gui_profiler_env_flag_accepts_common_truthy_values() {
        for value in ["1", "true", "TRUE", "yes", "on", " On "] {
            assert!(env_flag_value_enabled(Some(value)), "{value:?}");
        }
    }

    #[test]
    fn gui_profiler_env_flag_rejects_empty_or_false_values() {
        for value in [None, Some(""), Some("0"), Some("false"), Some("off")] {
            assert!(!env_flag_value_enabled(value), "{value:?}");
        }
    }
}
