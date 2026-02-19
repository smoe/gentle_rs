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
mod macos_native_help_menu {
    use std::{
        cell::RefCell,
        sync::atomic::{AtomicBool, Ordering},
    };

    use crate::app;
    use objc2::rc::Retained;
    use objc2::runtime::AnyObject;
    use objc2::{declare_class, msg_send_id, mutability, sel, ClassType, DeclaredClass};
    use objc2_app_kit::NSApplication;
    use objc2_foundation::{ns_string, MainThreadMarker, NSObject, NSObjectProtocol};

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
