//! About/help metadata and version presentation helpers.

pub const GENTLE_DISPLAY_VERSION: &str = env!("GENTLE_DISPLAY_VERSION");
pub const GENTLE_BUILD_N: &str = env!("GENTLE_BUILD_N");
pub const GENTLE_PACKAGE_VERSION: &str = env!("CARGO_PKG_VERSION");
pub const GENTLE_GIT_COMMIT: &str = env!("GENTLE_GIT_COMMIT");
pub const GENTLE_SOURCE_REVISION: &str = env!("GENTLE_SOURCE_REVISION");
pub const GENTLE_DESCRIPTION: &str = env!("CARGO_PKG_DESCRIPTION");
pub const GENTLE_REPOSITORY_URL: &str = env!("CARGO_PKG_REPOSITORY");
pub const GENTLE_DOCUMENTATION_URL: &str = "https://github.com/smoe/gentle_rs/tree/main/docs";
pub const GENTLE_ACKNOWLEDGEMENTS_URL: &str =
    "https://github.com/smoe/gentle_rs/blob/main/ACKNOWLEDGEMENTS.md";
pub const GENTLE_LICENSE: &str = "GPL-2.0-or-later";

pub fn version_cli_text() -> String {
    format!(
        "GENtle {}\nBuild {}\nSource revision {}\nCross-platform DNA cloning workbench",
        GENTLE_PACKAGE_VERSION, GENTLE_BUILD_N, GENTLE_SOURCE_REVISION
    )
}

pub fn about_clipboard_text() -> String {
    format!(
        "{}\n{}\nDisplay version {}\nRepository {}\nLicense {}",
        version_cli_text(),
        GENTLE_DESCRIPTION,
        GENTLE_DISPLAY_VERSION,
        GENTLE_REPOSITORY_URL,
        GENTLE_LICENSE
    )
}

#[cfg(all(feature = "desktop-gui", target_os = "macos"))]
pub fn install_native_about_menu_bridge() {
    macos_native_about_menu::install();
}

#[cfg(not(all(feature = "desktop-gui", target_os = "macos")))]
pub fn install_native_about_menu_bridge() {}

#[cfg(all(feature = "desktop-gui", target_os = "macos"))]
pub fn install_native_help_menu_bridge() {
    macos_native_help_menu::install();
}

#[cfg(not(all(feature = "desktop-gui", target_os = "macos")))]
pub fn install_native_help_menu_bridge() {}

#[cfg(all(feature = "desktop-gui", target_os = "macos"))]
pub fn install_native_settings_menu_bridge() {
    macos_native_settings_menu::install();
}

#[cfg(not(all(feature = "desktop-gui", target_os = "macos")))]
pub fn install_native_settings_menu_bridge() {}

#[cfg(all(feature = "desktop-gui", target_os = "macos"))]
pub fn install_native_windows_menu_bridge() {
    macos_native_windows_menu::install();
}

#[cfg(not(all(feature = "desktop-gui", target_os = "macos")))]
pub fn install_native_windows_menu_bridge() {}

#[cfg(all(feature = "desktop-gui", target_os = "macos"))]
pub fn install_native_app_windows_menu_bridge() {
    macos_native_app_windows_menu::install();
}

#[cfg(not(all(feature = "desktop-gui", target_os = "macos")))]
pub fn install_native_app_windows_menu_bridge() {}

#[cfg(all(feature = "desktop-gui", target_os = "macos"))]
pub fn sync_native_open_windows_menu(entries: &[(u64, String)], active_key: Option<u64>) {
    macos_native_windows_menu::sync_entries(entries, active_key);
    macos_native_app_windows_menu::sync_entries(entries, active_key);
}

#[cfg(not(all(feature = "desktop-gui", target_os = "macos")))]
pub fn sync_native_open_windows_menu(_entries: &[(u64, String)], _active_key: Option<u64>) {}

#[cfg(all(feature = "desktop-gui", target_os = "macos"))]
mod macos_native_about_menu {
    use std::{
        cell::RefCell,
        sync::atomic::{AtomicBool, Ordering},
    };

    use crate::{app, i18n};
    use objc2::rc::Retained;
    use objc2::runtime::AnyObject;
    use objc2::{MainThreadMarker, MainThreadOnly, define_class, msg_send, sel};
    use objc2_app_kit::NSApplication;
    use objc2_foundation::{NSObject, NSObjectProtocol, NSString};

    static INSTALLED: AtomicBool = AtomicBool::new(false);

    thread_local! {
        static ABOUT_MENU_TARGET: RefCell<Option<Retained<AboutMenuTarget>>> = const { RefCell::new(None) };
    }

    define_class!(
        #[unsafe(super(NSObject))]
        #[thread_kind = MainThreadOnly]
        #[name = "GENtleAboutMenuTarget"]
        #[ivars = ()]
        struct AboutMenuTarget;

        impl AboutMenuTarget {
            #[unsafe(method(openGentleAbout:))]
            fn open_gentle_about(&self, _sender: Option<&AnyObject>) {
                app::request_open_about_from_native_menu();
            }
        }

        unsafe impl NSObjectProtocol for AboutMenuTarget {}
    );

    impl AboutMenuTarget {
        fn new(mtm: MainThreadMarker) -> Retained<Self> {
            let this = Self::alloc(mtm).set_ivars(());
            unsafe { msg_send![super(this), init] }
        }
    }

    pub(super) fn install() {
        let Some(mtm) = MainThreadMarker::new() else {
            return;
        };
        let app = NSApplication::sharedApplication(mtm);
        let Some(main_menu) = app.mainMenu() else {
            return;
        };
        let Some(app_menu_item) = main_menu.itemAtIndex(0) else {
            return;
        };
        let Some(app_menu) = app_menu_item.submenu() else {
            return;
        };
        // Winit creates the standard About item first. Retarget that item so
        // macOS and GENtle's Help menu open the same application-owned panel.
        let Some(item) = app_menu.itemAtIndex(0) else {
            return;
        };
        item.setTitle(&NSString::from_str(&i18n::tr("menu.help.about")));

        if INSTALLED.load(Ordering::Relaxed) {
            return;
        }

        let target = AboutMenuTarget::new(mtm);
        let target_obj: &AnyObject = target.as_ref();
        unsafe {
            item.setTarget(Some(target_obj));
            item.setAction(Some(sel!(openGentleAbout:)));
        }
        ABOUT_MENU_TARGET.with(|slot| {
            *slot.borrow_mut() = Some(target);
        });
        INSTALLED.store(true, Ordering::Relaxed);
    }
}

#[cfg(all(feature = "desktop-gui", target_os = "macos"))]
mod macos_native_help_menu {
    use std::{
        cell::RefCell,
        sync::atomic::{AtomicBool, Ordering},
    };

    use crate::app;
    use objc2::rc::Retained;
    use objc2::runtime::AnyObject;
    use objc2::{MainThreadMarker, MainThreadOnly, define_class, msg_send, sel};
    use objc2_app_kit::NSApplication;
    use objc2_foundation::{NSObject, NSObjectProtocol, ns_string};

    static INSTALLED: AtomicBool = AtomicBool::new(false);

    thread_local! {
        static HELP_MENU_TARGET: RefCell<Option<Retained<HelpMenuTarget>>> = const { RefCell::new(None) };
    }

    define_class!(
        #[unsafe(super(NSObject))]
        #[thread_kind = MainThreadOnly]
        #[name = "GENtleHelpMenuTarget"]
        #[ivars = ()]
        struct HelpMenuTarget;

        impl HelpMenuTarget {
            #[unsafe(method(openGentleHelp:))]
            fn open_gentle_help(&self, _sender: Option<&AnyObject>) {
                app::request_open_help_from_native_menu();
            }
        }

        unsafe impl NSObjectProtocol for HelpMenuTarget {}
    );

    impl HelpMenuTarget {
        fn new(mtm: MainThreadMarker) -> Retained<Self> {
            let this = Self::alloc(mtm).set_ivars(());
            unsafe { msg_send![super(this), init] }
        }
    }

    pub(super) fn install() {
        if INSTALLED.load(Ordering::Relaxed) {
            return;
        }

        let Some(mtm) = MainThreadMarker::new() else {
            return;
        };
        let app = NSApplication::sharedApplication(mtm);
        let Some(main_menu) = app.mainMenu() else {
            return;
        };
        let Some(app_menu_item) = main_menu.itemAtIndex(0) else {
            return;
        };
        let Some(app_menu) = app_menu_item.submenu() else {
            return;
        };

        let item_title = ns_string!("GENtle Help...");
        if app_menu.indexOfItemWithTitle(item_title) >= 0 {
            INSTALLED.store(true, Ordering::Relaxed);
            return;
        }

        let item = unsafe {
            app_menu.insertItemWithTitle_action_keyEquivalent_atIndex(
                item_title,
                Some(sel!(openGentleHelp:)),
                ns_string!(""),
                1,
            )
        };
        let target = HelpMenuTarget::new(mtm);
        let target_obj: &AnyObject = target.as_ref();
        unsafe {
            item.setTarget(Some(target_obj));
        }
        HELP_MENU_TARGET.with(|slot| {
            *slot.borrow_mut() = Some(target);
        });
        INSTALLED.store(true, Ordering::Relaxed);
    }
}

#[cfg(all(feature = "desktop-gui", target_os = "macos"))]
mod macos_native_settings_menu {
    use std::{
        cell::RefCell,
        sync::atomic::{AtomicBool, Ordering},
    };

    use crate::app;
    use objc2::rc::Retained;
    use objc2::runtime::AnyObject;
    use objc2::{MainThreadMarker, MainThreadOnly, define_class, msg_send, sel};
    use objc2_app_kit::NSApplication;
    use objc2_foundation::{NSObject, NSObjectProtocol, ns_string};

    static INSTALLED: AtomicBool = AtomicBool::new(false);

    thread_local! {
        static SETTINGS_MENU_TARGET: RefCell<Option<Retained<SettingsMenuTarget>>> = const { RefCell::new(None) };
    }

    define_class!(
        #[unsafe(super(NSObject))]
        #[thread_kind = MainThreadOnly]
        #[name = "GENtleSettingsMenuTarget"]
        #[ivars = ()]
        struct SettingsMenuTarget;

        impl SettingsMenuTarget {
            #[unsafe(method(openGentleSettings:))]
            fn open_gentle_settings(&self, _sender: Option<&AnyObject>) {
                app::request_open_settings_from_native_menu();
            }
        }

        unsafe impl NSObjectProtocol for SettingsMenuTarget {}
    );

    impl SettingsMenuTarget {
        fn new(mtm: MainThreadMarker) -> Retained<Self> {
            let this = Self::alloc(mtm).set_ivars(());
            unsafe { msg_send![super(this), init] }
        }
    }

    pub(super) fn install() {
        if INSTALLED.load(Ordering::Relaxed) {
            return;
        }

        let Some(mtm) = MainThreadMarker::new() else {
            return;
        };
        let app = NSApplication::sharedApplication(mtm);
        let Some(main_menu) = app.mainMenu() else {
            return;
        };
        let Some(app_menu_item) = main_menu.itemAtIndex(0) else {
            return;
        };
        let Some(app_menu) = app_menu_item.submenu() else {
            return;
        };

        let item_title = ns_string!("GENtle Settings...");
        if app_menu.indexOfItemWithTitle(item_title) >= 0 {
            INSTALLED.store(true, Ordering::Relaxed);
            return;
        }

        let item = unsafe {
            app_menu.insertItemWithTitle_action_keyEquivalent_atIndex(
                item_title,
                Some(sel!(openGentleSettings:)),
                ns_string!(","),
                2,
            )
        };
        let target = SettingsMenuTarget::new(mtm);
        let target_obj: &AnyObject = target.as_ref();
        unsafe {
            item.setTarget(Some(target_obj));
        }
        SETTINGS_MENU_TARGET.with(|slot| {
            *slot.borrow_mut() = Some(target);
        });
        INSTALLED.store(true, Ordering::Relaxed);
    }
}

#[cfg(all(feature = "desktop-gui", target_os = "macos"))]
mod macos_native_windows_menu {
    use std::{
        cell::RefCell,
        sync::atomic::{AtomicBool, Ordering},
    };

    use crate::app;
    use objc2::rc::Retained;
    use objc2::runtime::AnyObject;
    use objc2::{MainThreadMarker, MainThreadOnly, define_class, msg_send, sel};
    use objc2_app_kit::{NSApplication, NSMenu};
    use objc2_foundation::{NSObject, NSObjectProtocol, NSString, ns_string};

    static INSTALLED: AtomicBool = AtomicBool::new(false);

    thread_local! {
        static WINDOWS_MENU_TARGET: RefCell<Option<Retained<WindowsMenuTarget>>> = const { RefCell::new(None) };
    }

    define_class!(
        #[unsafe(super(NSObject))]
        #[thread_kind = MainThreadOnly]
        #[name = "GENtleWindowsMenuTarget"]
        #[ivars = ()]
        struct WindowsMenuTarget;

        impl WindowsMenuTarget {
            #[unsafe(method(openGentleWindows:))]
            fn open_gentle_windows(&self, _sender: Option<&AnyObject>) {
                app::request_open_windows_from_native_menu();
            }

            #[unsafe(method(focusGentleWindowEntry:))]
            fn focus_gentle_window_entry(&self, sender: Option<&AnyObject>) {
                let Some(sender) = sender else {
                    app::request_open_windows_from_native_menu();
                    return;
                };
                let raw_tag: isize = unsafe { msg_send![sender, tag] };
                if raw_tag >= 0 {
                    app::request_focus_window_key_from_native_menu(raw_tag as u64);
                } else {
                    app::request_open_windows_from_native_menu();
                }
            }
        }

        unsafe impl NSObjectProtocol for WindowsMenuTarget {}
    );

    impl WindowsMenuTarget {
        fn new(mtm: MainThreadMarker) -> Retained<Self> {
            let this = Self::alloc(mtm).set_ivars(());
            unsafe { msg_send![super(this), init] }
        }
    }

    fn find_windows_item(window_menu: &NSMenu) -> Option<Retained<objc2_app_kit::NSMenuItem>> {
        let item_title = ns_string!("GENtle Open Windows…");
        let index = window_menu.indexOfItemWithTitle(item_title);
        if index < 0 {
            return None;
        }
        window_menu.itemAtIndex(index)
    }

    fn try_update_existing_entry_states(
        submenu: &NSMenu,
        entries: &[(u64, String)],
        active_key: Option<u64>,
    ) -> bool {
        if entries.is_empty() {
            return false;
        }
        let item_count = submenu.numberOfItems();
        if item_count < 0 || item_count as usize != entries.len() {
            return false;
        }
        for (idx, (key, _title)) in entries.iter().enumerate() {
            let Some(item) = submenu.itemAtIndex(idx as isize) else {
                return false;
            };
            let existing_tag: isize = unsafe { msg_send![&item, tag] };
            if existing_tag != *key as isize {
                return false;
            }
        }
        for (idx, (key, _title)) in entries.iter().enumerate() {
            let Some(item) = submenu.itemAtIndex(idx as isize) else {
                return false;
            };
            let state_value: isize = if Some(*key) == active_key { 1 } else { 0 };
            let _: () = unsafe { msg_send![&item, setState: state_value] };
            item.setEnabled(true);
        }
        true
    }

    pub(super) fn install() {
        if INSTALLED.load(Ordering::Relaxed) {
            return;
        }

        let Some(mtm) = MainThreadMarker::new() else {
            return;
        };
        let app = NSApplication::sharedApplication(mtm);
        let Some(main_menu) = app.mainMenu() else {
            return;
        };
        let Some(window_menu_item) = main_menu.itemWithTitle(ns_string!("Window")) else {
            return;
        };
        let Some(window_menu) = window_menu_item.submenu() else {
            return;
        };

        let item_title = ns_string!("GENtle Open Windows…");
        if window_menu.indexOfItemWithTitle(item_title) >= 0 {
            INSTALLED.store(true, Ordering::Relaxed);
            return;
        }

        let item = unsafe {
            window_menu.insertItemWithTitle_action_keyEquivalent_atIndex(
                item_title,
                Some(sel!(openGentleWindows:)),
                ns_string!("`"),
                0,
            )
        };
        let target = WindowsMenuTarget::new(mtm);
        let target_obj: &AnyObject = target.as_ref();
        unsafe {
            item.setTarget(Some(target_obj));
        }
        WINDOWS_MENU_TARGET.with(|slot| {
            *slot.borrow_mut() = Some(target);
        });
        INSTALLED.store(true, Ordering::Relaxed);
    }

    pub(super) fn sync_entries(entries: &[(u64, String)], active_key: Option<u64>) {
        let Some(mtm) = MainThreadMarker::new() else {
            return;
        };
        let app = NSApplication::sharedApplication(mtm);
        let Some(main_menu) = app.mainMenu() else {
            return;
        };
        let Some(window_menu_item) = main_menu.itemWithTitle(ns_string!("Window")) else {
            return;
        };
        let Some(window_menu) = window_menu_item.submenu() else {
            return;
        };
        let Some(item) = find_windows_item(&window_menu) else {
            return;
        };
        let submenu = if let Some(existing) = item.submenu() {
            existing
        } else {
            let created = NSMenu::initWithTitle(mtm.alloc(), ns_string!("GENtle Open Windows"));
            item.setSubmenu(Some(&created));
            created
        };
        if try_update_existing_entry_states(&submenu, entries, active_key) {
            return;
        }
        submenu.removeAllItems();
        if entries.is_empty() {
            let entry = unsafe {
                submenu.addItemWithTitle_action_keyEquivalent(
                    ns_string!("No open windows"),
                    None,
                    ns_string!(""),
                )
            };
            entry.setEnabled(false);
            return;
        }

        for (key, title) in entries {
            let title = NSString::from_str(title);
            let entry = unsafe {
                submenu.addItemWithTitle_action_keyEquivalent(&title, None, ns_string!(""))
            };
            WINDOWS_MENU_TARGET.with(|slot| {
                if let Some(target) = slot.borrow().as_ref() {
                    let target_obj: &AnyObject = target.as_ref();
                    unsafe {
                        entry.setTarget(Some(target_obj));
                        entry.setAction(Some(sel!(focusGentleWindowEntry:)));
                    }
                    entry.setTag(*key as isize);
                    let state_value: isize = if Some(*key) == active_key { 1 } else { 0 };
                    let _: () = unsafe { msg_send![&entry, setState: state_value] };
                    entry.setEnabled(true);
                } else {
                    entry.setEnabled(false);
                }
            });
        }
    }
}

#[cfg(all(feature = "desktop-gui", target_os = "macos"))]
mod macos_native_app_windows_menu {
    use std::{
        cell::RefCell,
        sync::atomic::{AtomicBool, Ordering},
    };

    use crate::app;
    use objc2::rc::Retained;
    use objc2::runtime::AnyObject;
    use objc2::{MainThreadMarker, MainThreadOnly, define_class, msg_send, sel};
    use objc2_app_kit::{NSApplication, NSMenu};
    use objc2_foundation::{NSObject, NSObjectProtocol, NSString, ns_string};

    static INSTALLED: AtomicBool = AtomicBool::new(false);

    thread_local! {
        static APP_WINDOWS_MENU_TARGET: RefCell<Option<Retained<AppWindowsMenuTarget>>> = const { RefCell::new(None) };
    }

    define_class!(
        #[unsafe(super(NSObject))]
        #[thread_kind = MainThreadOnly]
        #[name = "GENtleAppWindowsMenuTarget"]
        #[ivars = ()]
        struct AppWindowsMenuTarget;

        impl AppWindowsMenuTarget {
            #[unsafe(method(openGentleWindowsFromAppMenu:))]
            fn open_gentle_windows_from_app_menu(&self, _sender: Option<&AnyObject>) {
                app::request_open_windows_from_native_menu();
            }

            #[unsafe(method(focusGentleWindowEntryFromAppMenu:))]
            fn focus_gentle_window_entry_from_app_menu(&self, sender: Option<&AnyObject>) {
                let Some(sender) = sender else {
                    app::request_open_windows_from_native_menu();
                    return;
                };
                let raw_tag: isize = unsafe { msg_send![sender, tag] };
                if raw_tag >= 0 {
                    app::request_focus_window_key_from_native_menu(raw_tag as u64);
                } else {
                    app::request_open_windows_from_native_menu();
                }
            }
        }

        unsafe impl NSObjectProtocol for AppWindowsMenuTarget {}
    );

    impl AppWindowsMenuTarget {
        fn new(mtm: MainThreadMarker) -> Retained<Self> {
            let this = Self::alloc(mtm).set_ivars(());
            unsafe { msg_send![super(this), init] }
        }
    }

    fn find_app_windows_item(app_menu: &NSMenu) -> Option<Retained<objc2_app_kit::NSMenuItem>> {
        let item_title = ns_string!("GENtle Windows…");
        let index = app_menu.indexOfItemWithTitle(item_title);
        if index < 0 {
            return None;
        }
        app_menu.itemAtIndex(index)
    }

    fn try_update_existing_entry_states(
        submenu: &NSMenu,
        entries: &[(u64, String)],
        active_key: Option<u64>,
    ) -> bool {
        if entries.is_empty() {
            return false;
        }
        let item_count = submenu.numberOfItems();
        if item_count < 0 || item_count as usize != entries.len() {
            return false;
        }
        for (idx, (key, _title)) in entries.iter().enumerate() {
            let Some(item) = submenu.itemAtIndex(idx as isize) else {
                return false;
            };
            let existing_tag: isize = unsafe { msg_send![&item, tag] };
            if existing_tag != *key as isize {
                return false;
            }
        }
        for (idx, (key, _title)) in entries.iter().enumerate() {
            let Some(item) = submenu.itemAtIndex(idx as isize) else {
                return false;
            };
            let state_value: isize = if Some(*key) == active_key { 1 } else { 0 };
            let _: () = unsafe { msg_send![&item, setState: state_value] };
            item.setEnabled(true);
        }
        true
    }

    pub(super) fn install() {
        if INSTALLED.load(Ordering::Relaxed) {
            return;
        }

        let Some(mtm) = MainThreadMarker::new() else {
            return;
        };
        let app = NSApplication::sharedApplication(mtm);
        let Some(main_menu) = app.mainMenu() else {
            return;
        };
        let Some(app_menu_item) = main_menu.itemAtIndex(0) else {
            return;
        };
        let Some(app_menu) = app_menu_item.submenu() else {
            return;
        };

        let item_title = ns_string!("GENtle Windows…");
        if app_menu.indexOfItemWithTitle(item_title) >= 0 {
            INSTALLED.store(true, Ordering::Relaxed);
            return;
        }

        let item = unsafe {
            app_menu.insertItemWithTitle_action_keyEquivalent_atIndex(
                item_title,
                Some(sel!(openGentleWindowsFromAppMenu:)),
                ns_string!(""),
                3,
            )
        };
        let target = AppWindowsMenuTarget::new(mtm);
        let target_obj: &AnyObject = target.as_ref();
        unsafe {
            item.setTarget(Some(target_obj));
        }
        APP_WINDOWS_MENU_TARGET.with(|slot| {
            *slot.borrow_mut() = Some(target);
        });
        INSTALLED.store(true, Ordering::Relaxed);
    }

    pub(super) fn sync_entries(entries: &[(u64, String)], active_key: Option<u64>) {
        let Some(mtm) = MainThreadMarker::new() else {
            return;
        };
        let app = NSApplication::sharedApplication(mtm);
        let Some(main_menu) = app.mainMenu() else {
            return;
        };
        let Some(app_menu_item) = main_menu.itemAtIndex(0) else {
            return;
        };
        let Some(app_menu) = app_menu_item.submenu() else {
            return;
        };
        let Some(item) = find_app_windows_item(&app_menu) else {
            return;
        };
        let submenu = if let Some(existing) = item.submenu() {
            existing
        } else {
            let created = NSMenu::initWithTitle(mtm.alloc(), ns_string!("GENtle Windows"));
            item.setSubmenu(Some(&created));
            created
        };
        if try_update_existing_entry_states(&submenu, entries, active_key) {
            return;
        }
        submenu.removeAllItems();
        if entries.is_empty() {
            let entry = unsafe {
                submenu.addItemWithTitle_action_keyEquivalent(
                    ns_string!("No open windows"),
                    None,
                    ns_string!(""),
                )
            };
            entry.setEnabled(false);
            return;
        }

        for (key, title) in entries {
            let title = NSString::from_str(title);
            let entry = unsafe {
                submenu.addItemWithTitle_action_keyEquivalent(&title, None, ns_string!(""))
            };
            APP_WINDOWS_MENU_TARGET.with(|slot| {
                if let Some(target) = slot.borrow().as_ref() {
                    let target_obj: &AnyObject = target.as_ref();
                    unsafe {
                        entry.setTarget(Some(target_obj));
                        entry.setAction(Some(sel!(focusGentleWindowEntryFromAppMenu:)));
                    }
                    entry.setTag(*key as isize);
                    let state_value: isize = if Some(*key) == active_key { 1 } else { 0 };
                    let _: () = unsafe { msg_send![&entry, setState: state_value] };
                    entry.setEnabled(true);
                } else {
                    entry.setEnabled(false);
                }
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn about_metadata_is_bound_to_the_built_package() {
        assert!(!GENTLE_PACKAGE_VERSION.is_empty());
        assert!(!GENTLE_BUILD_N.is_empty());
        assert!(!GENTLE_SOURCE_REVISION.is_empty());
        if GENTLE_GIT_COMMIT == "unknown" {
            assert_eq!(GENTLE_SOURCE_REVISION, GENTLE_PACKAGE_VERSION);
        } else {
            assert!((7..=64).contains(&GENTLE_GIT_COMMIT.len()));
            assert!(
                GENTLE_GIT_COMMIT
                    .bytes()
                    .all(|byte| byte.is_ascii_hexdigit())
            );
            assert_eq!(
                GENTLE_SOURCE_REVISION,
                format!("{GENTLE_PACKAGE_VERSION}+git.{GENTLE_GIT_COMMIT}")
            );
        }
        assert!(GENTLE_REPOSITORY_URL.starts_with("https://github.com/"));
        assert_eq!(GENTLE_LICENSE, "GPL-2.0-or-later");

        let clipboard = about_clipboard_text();
        assert!(clipboard.contains(GENTLE_PACKAGE_VERSION));
        assert!(clipboard.contains(GENTLE_BUILD_N));
        assert!(clipboard.contains(GENTLE_SOURCE_REVISION));
        assert!(clipboard.contains(GENTLE_REPOSITORY_URL));
        assert!(clipboard.contains(GENTLE_LICENSE));
    }
}
