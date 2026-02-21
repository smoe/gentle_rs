pub const GENTLE_DISPLAY_VERSION: &str = env!("GENTLE_DISPLAY_VERSION");
pub const GENTLE_BUILD_N: &str = env!("GENTLE_BUILD_N");

pub fn version_cli_text() -> String {
    format!(
        "GENtle {}\nBuild {}\nCross-platform DNA cloning workbench",
        GENTLE_DISPLAY_VERSION, GENTLE_BUILD_N
    )
}

#[cfg(target_os = "macos")]
pub fn show_native_about_panel() -> bool {
    use objc2_app_kit::NSApplication;
    use objc2_foundation::MainThreadMarker;

    let Some(mtm) = MainThreadMarker::new() else {
        return false;
    };
    let app = NSApplication::sharedApplication(mtm);
    unsafe {
        app.orderFrontStandardAboutPanel(None);
    }
    true
}

#[cfg(not(target_os = "macos"))]
pub fn show_native_about_panel() -> bool {
    false
}

#[cfg(target_os = "macos")]
pub fn install_native_help_menu_bridge() {
    macos_native_help_menu::install();
}

#[cfg(not(target_os = "macos"))]
pub fn install_native_help_menu_bridge() {}

#[cfg(target_os = "macos")]
pub fn install_native_settings_menu_bridge() {
    macos_native_settings_menu::install();
}

#[cfg(not(target_os = "macos"))]
pub fn install_native_settings_menu_bridge() {}

#[cfg(target_os = "macos")]
pub fn install_native_windows_menu_bridge() {
    macos_native_windows_menu::install();
}

#[cfg(not(target_os = "macos"))]
pub fn install_native_windows_menu_bridge() {}

#[cfg(target_os = "macos")]
pub fn install_native_app_windows_menu_bridge() {
    macos_native_app_windows_menu::install();
}

#[cfg(not(target_os = "macos"))]
pub fn install_native_app_windows_menu_bridge() {}

#[cfg(target_os = "macos")]
mod macos_native_help_menu {
    use std::{
        cell::RefCell,
        sync::atomic::{AtomicBool, Ordering},
    };

    use crate::app;
    use objc2::rc::Retained;
    use objc2::runtime::AnyObject;
    use objc2::{ClassType, DeclaredClass, declare_class, msg_send_id, mutability, sel};
    use objc2_app_kit::NSApplication;
    use objc2_foundation::{MainThreadMarker, NSObject, NSObjectProtocol, ns_string};

    static INSTALLED: AtomicBool = AtomicBool::new(false);

    thread_local! {
        static HELP_MENU_TARGET: RefCell<Option<Retained<HelpMenuTarget>>> = const { RefCell::new(None) };
    }

    declare_class!(
        struct HelpMenuTarget;

        unsafe impl ClassType for HelpMenuTarget {
            type Super = NSObject;
            type Mutability = mutability::MainThreadOnly;
            const NAME: &'static str = "GENtleHelpMenuTarget";
        }

        impl DeclaredClass for HelpMenuTarget {
            type Ivars = ();
        }

        unsafe impl NSObjectProtocol for HelpMenuTarget {}
        unsafe impl HelpMenuTarget {
            #[method(openGentleHelp:)]
            fn open_gentle_help(&self, _sender: Option<&AnyObject>) {
                app::request_open_help_from_native_menu();
            }
        }
    );

    impl HelpMenuTarget {
        fn new(mtm: MainThreadMarker) -> Retained<Self> {
            let this = mtm.alloc();
            let this = this.set_ivars(());
            unsafe { msg_send_id![super(this), init] }
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
        let Some(main_menu) = (unsafe { app.mainMenu() }) else {
            return;
        };
        let Some(app_menu_item) = (unsafe { main_menu.itemAtIndex(0) }) else {
            return;
        };
        let Some(app_menu) = (unsafe { app_menu_item.submenu() }) else {
            return;
        };

        let item_title = ns_string!("GENtle Help...");
        if unsafe { app_menu.indexOfItemWithTitle(item_title) } >= 0 {
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

#[cfg(target_os = "macos")]
mod macos_native_settings_menu {
    use std::{
        cell::RefCell,
        sync::atomic::{AtomicBool, Ordering},
    };

    use crate::app;
    use objc2::rc::Retained;
    use objc2::runtime::AnyObject;
    use objc2::{ClassType, DeclaredClass, declare_class, msg_send_id, mutability, sel};
    use objc2_app_kit::NSApplication;
    use objc2_foundation::{MainThreadMarker, NSObject, NSObjectProtocol, ns_string};

    static INSTALLED: AtomicBool = AtomicBool::new(false);

    thread_local! {
        static SETTINGS_MENU_TARGET: RefCell<Option<Retained<SettingsMenuTarget>>> = const { RefCell::new(None) };
    }

    declare_class!(
        struct SettingsMenuTarget;

        unsafe impl ClassType for SettingsMenuTarget {
            type Super = NSObject;
            type Mutability = mutability::MainThreadOnly;
            const NAME: &'static str = "GENtleSettingsMenuTarget";
        }

        impl DeclaredClass for SettingsMenuTarget {
            type Ivars = ();
        }

        unsafe impl NSObjectProtocol for SettingsMenuTarget {}
        unsafe impl SettingsMenuTarget {
            #[method(openGentleSettings:)]
            fn open_gentle_settings(&self, _sender: Option<&AnyObject>) {
                app::request_open_settings_from_native_menu();
            }
        }
    );

    impl SettingsMenuTarget {
        fn new(mtm: MainThreadMarker) -> Retained<Self> {
            let this = mtm.alloc();
            let this = this.set_ivars(());
            unsafe { msg_send_id![super(this), init] }
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
        let Some(main_menu) = (unsafe { app.mainMenu() }) else {
            return;
        };
        let Some(app_menu_item) = (unsafe { main_menu.itemAtIndex(0) }) else {
            return;
        };
        let Some(app_menu) = (unsafe { app_menu_item.submenu() }) else {
            return;
        };

        let item_title = ns_string!("GENtle Settings...");
        if unsafe { app_menu.indexOfItemWithTitle(item_title) } >= 0 {
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

#[cfg(target_os = "macos")]
mod macos_native_windows_menu {
    use std::{
        cell::RefCell,
        sync::atomic::{AtomicBool, Ordering},
    };

    use crate::app;
    use objc2::rc::Retained;
    use objc2::runtime::AnyObject;
    use objc2::{ClassType, DeclaredClass, declare_class, msg_send_id, mutability, sel};
    use objc2_app_kit::NSApplication;
    use objc2_foundation::{MainThreadMarker, NSObject, NSObjectProtocol, ns_string};

    static INSTALLED: AtomicBool = AtomicBool::new(false);

    thread_local! {
        static WINDOWS_MENU_TARGET: RefCell<Option<Retained<WindowsMenuTarget>>> = const { RefCell::new(None) };
    }

    declare_class!(
        struct WindowsMenuTarget;

        unsafe impl ClassType for WindowsMenuTarget {
            type Super = NSObject;
            type Mutability = mutability::MainThreadOnly;
            const NAME: &'static str = "GENtleWindowsMenuTarget";
        }

        impl DeclaredClass for WindowsMenuTarget {
            type Ivars = ();
        }

        unsafe impl NSObjectProtocol for WindowsMenuTarget {}
        unsafe impl WindowsMenuTarget {
            #[method(openGentleWindows:)]
            fn open_gentle_windows(&self, _sender: Option<&AnyObject>) {
                app::request_open_windows_from_native_menu();
            }
        }
    );

    impl WindowsMenuTarget {
        fn new(mtm: MainThreadMarker) -> Retained<Self> {
            let this = mtm.alloc();
            let this = this.set_ivars(());
            unsafe { msg_send_id![super(this), init] }
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
        let Some(main_menu) = (unsafe { app.mainMenu() }) else {
            return;
        };
        let Some(window_menu_item) = (unsafe { main_menu.itemWithTitle(ns_string!("Window")) })
        else {
            return;
        };
        let Some(window_menu) = (unsafe { window_menu_item.submenu() }) else {
            return;
        };

        let item_title = ns_string!("GENtle Open Windows…");
        if unsafe { window_menu.indexOfItemWithTitle(item_title) } >= 0 {
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
}

#[cfg(target_os = "macos")]
mod macos_native_app_windows_menu {
    use std::{
        cell::RefCell,
        sync::atomic::{AtomicBool, Ordering},
    };

    use crate::app;
    use objc2::rc::Retained;
    use objc2::runtime::AnyObject;
    use objc2::{ClassType, DeclaredClass, declare_class, msg_send_id, mutability, sel};
    use objc2_app_kit::NSApplication;
    use objc2_foundation::{MainThreadMarker, NSObject, NSObjectProtocol, ns_string};

    static INSTALLED: AtomicBool = AtomicBool::new(false);

    thread_local! {
        static APP_WINDOWS_MENU_TARGET: RefCell<Option<Retained<AppWindowsMenuTarget>>> = const { RefCell::new(None) };
    }

    declare_class!(
        struct AppWindowsMenuTarget;

        unsafe impl ClassType for AppWindowsMenuTarget {
            type Super = NSObject;
            type Mutability = mutability::MainThreadOnly;
            const NAME: &'static str = "GENtleAppWindowsMenuTarget";
        }

        impl DeclaredClass for AppWindowsMenuTarget {
            type Ivars = ();
        }

        unsafe impl NSObjectProtocol for AppWindowsMenuTarget {}
        unsafe impl AppWindowsMenuTarget {
            #[method(openGentleWindowsFromAppMenu:)]
            fn open_gentle_windows_from_app_menu(&self, _sender: Option<&AnyObject>) {
                app::request_open_windows_from_native_menu();
            }
        }
    );

    impl AppWindowsMenuTarget {
        fn new(mtm: MainThreadMarker) -> Retained<Self> {
            let this = mtm.alloc();
            let this = this.set_ivars(());
            unsafe { msg_send_id![super(this), init] }
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
        let Some(main_menu) = (unsafe { app.mainMenu() }) else {
            return;
        };
        let Some(app_menu_item) = (unsafe { main_menu.itemAtIndex(0) }) else {
            return;
        };
        let Some(app_menu) = (unsafe { app_menu_item.submenu() }) else {
            return;
        };

        let item_title = ns_string!("GENtle Windows…");
        if unsafe { app_menu.indexOfItemWithTitle(item_title) } >= 0 {
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
}
